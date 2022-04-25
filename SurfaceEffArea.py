#!/usr/bin/env python

# Written by: Alan Coleman
# Calculates the geometric aperture for surface/in-ice coincident events

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from icecube import dataclasses, icetray, dataio, radcube, rock_bottom
from icecube.icetray import I3Units

from sa_gen2_utils.PolygonTools import GetPolygonIntersection, GetHull, GetHullArea
from sa_gen2_utils.GeometryParser import GetDetectorGeometries, Get2DProjection, ProjectToObslev

import numpy as np
import os

ABS_PATH_HERE = str(os.path.dirname(os.path.realpath(__file__)))

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", type=str, help="Input data files.")
parser.add_argument("--output", type=str, default=ABS_PATH_HERE + "/../plots/SurfaceEffArea.pdf", help="Output data files.")
args = parser.parse_args()


scints, ants, doms = GetDetectorGeometries(args.input)

print("Read in", len(scints), "scints from the GCD file")
print("Read in", len(doms), "in-ice doms from the GCD file")

# make 3D by hand since this will fail otherwise
surfaceHull2D = GetHull(scints, 2)
surfaceHull3D = np.array([[val[0], val[1], scints[0][2]] for val in surfaceHull2D])
iniceHull = GetHull(doms, 3)

print("Surface array area:", GetHullArea(surfaceHull2D) / I3Units.km2, " km^2")

minX = min([pos.x for pos in scints])
maxX = max([pos.x for pos in scints])
minY = min([pos.y for pos in scints])
maxY = max([pos.y for pos in scints])

# Values to define the numerical integration
# The bins are a bit large, but this already takes a long time
dcos2 = 0.015
mincos2 = np.cos(70 * I3Units.degree) ** 2
maxcos2 = 1 - dcos2 / 2.0  # Do not remove the -dcos/2
dazi = 5 * I3Units.degree
dcore = 25 * I3Units.m

x_list = np.linspace(minX, maxX, int((maxX - minX) / dcore))
y_list = np.linspace(minY, maxY, int((maxY - minY) / dcore))
xs, ys = np.meshgrid(x_list, y_list)
xs = xs.flatten()
ys = ys.flatten()


effArea = np.zeros(shape=(len(xs)))
maxZen = np.zeros(shape=(len(xs)))

nAziBins = int(np.pi * 2 / dazi + 1)
aziValues = np.linspace(0.0, 2 * np.pi, nAziBins)

nCosBins = int((maxcos2 - mincos2) / dcos2 + 1)
cos2Values = np.linspace(mincos2, maxcos2, nCosBins)
effAreaZenith = np.zeros(shape=(nCosBins))

# These are the actual widths after the int process
dcos2 = cos2Values[1] - cos2Values[0]
dazi = aziValues[1] - aziValues[0]

effUnit = (x_list[1] - x_list[0]) * (y_list[1] - y_list[0]) * 0.5 * dcos2 * dazi

# Numerically intergrate the aperture

for icos2, cos2 in enumerate(cos2Values):
    zenith = np.arccos(np.sqrt(cos2))
    print("Calculating zenith angle: {0:0.2f}".format(zenith / I3Units.degree))
    for azi in aziValues:

        direction = dataclasses.I3Direction(zenith, azi)

        # Project the hulls into the shower CS and find overlap of two components
        inIce2DProj = Get2DProjection(direction, iniceHull)
        surf2DProj = Get2DProjection(direction, surfaceHull3D)
        intersection = GetPolygonIntersection(surf2DProj, inIce2DProj)

        # Skip if there isn't any intersection
        if not len(intersection):
            continue

        effAreaZenith[icos2] += GetHullArea(intersection)

        # Project the intersection back to the surface
        intersectionSurf = [ProjectToObslev(dataclasses.I3Position(pos[0], pos[1], 0), direction, scints[0][2]) for pos in intersection]
        intersectionSurf = [[pos[0], pos[1]] for pos in intersectionSurf]

        minx = min([pos[0] for pos in intersectionSurf])
        maxx = max([pos[0] for pos in intersectionSurf])
        miny = min([pos[1] for pos in intersectionSurf])
        maxy = max([pos[1] for pos in intersectionSurf])

        intersectionSurfPath = matplotlib.path.Path(intersectionSurf)

        for i, (x, y) in enumerate(zip(xs, ys)):

            if (not minx <= x <= maxx) or (not miny <= y <= maxy):
                continue
            if not intersectionSurfPath.contains_point((x, y)):
                continue

            if zenith == 0:
                effArea[i] += effUnit / 2.0  # Zenith step is only 0 to dcos/2
            else:
                effArea[i] += effUnit

            maxZen[i] = max([maxZen[i], zenith])


effAreaZenith /= nAziBins



NRows = 2
NCols = 2
gs = gridspec.GridSpec(NRows, NCols, wspace=0.3, hspace=0.3)
fig = plt.figure(figsize=(NCols * 6 * 1.2, NRows * 5 * 1.2))

### Plot of the geometric area of each pixel

ax = fig.add_subplot(gs[0])
ax.set_aspect("equal")
ax.set_xlabel("East / m")
ax.set_ylabel("North / m")

z = effArea.reshape(len(y_list), len(x_list))
z = np.flipud(z)
z = np.ma.masked_where(z == 0, z)
plot = ax.imshow(z / I3Units.m ** 2, extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)))
cbar = plt.colorbar(plot)
cbar.set_label(r"Pixel Aperture / m$^2$ sr")
ax.set_title(r"Geometric Aperture: {0:0.1f} km$^2$ sr".format(sum(effArea) / I3Units.km ** 2))


### Plot of the max zenith angle

ax = fig.add_subplot(gs[1])
ax.set_aspect("equal")

z = maxZen.reshape(len(y_list), len(x_list))
z = np.flipud(z)
z = np.ma.masked_where(z == 0, z)
plot = ax.imshow(z / I3Units.degree, extent=(np.amin(x_list), np.amax(x_list), np.amin(y_list), np.amax(y_list)))
ax.set_xlabel("East / m")
ax.set_ylabel("North / m")
cbar = plt.colorbar(plot)

minbar = min([val / I3Units.degree for val in maxZen if val > 0])
cbar.set_clim(minbar - 1, max(maxZen) / I3Units.degree + 1)
cbar.set_label("Maximum Zenith Angle / deg")

ax.set_title(" ")


### Declination dependence of the aperture
ax = fig.add_subplot(gs[2])

zenAngles = np.arccos(np.sqrt(cos2Values))
ax.plot(np.sin(-np.pi / 2.0 + zenAngles), effAreaZenith / I3Units.km ** 2)
ax.scatter(np.sin(-np.pi / 2.0 + zenAngles), effAreaZenith / I3Units.km ** 2)

ax.set_xlabel(r"sin($\delta$)")
ax.set_ylabel(r"Effective Area (km$^2$)")
ax.set_title("Declination Dependence")


ax = fig.add_subplot(gs[3])

ax.plot(np.cos(zenAngles), effAreaZenith / I3Units.km ** 2)
ax.scatter(np.cos(zenAngles), effAreaZenith / I3Units.km ** 2)

ax.set_xlabel(r"cos($\theta$)")
ax.set_ylabel(r"Effective Area (km$^2$)")
ax.set_title("Zenith Dependence")

print("Makeing file", args.output)
plt.savefig(args.output, bbox_inches="tight")
plt.close()
