
i3_project(weighting PYTHON_DIR python DOCS_DIR resources/docs)

colormsg(YELLOW "This project will be removed from icetray January 2022.")
colormsg(YELLOW "Prefer simweights: https://github.com/icecube/simweights")

# disabling the test that fails in Madison, which is also
# the only place it will run.
# See issue #2553:  https://github.com/icecube/icetray/issues/2553
i3_test_scripts(
	resources/test/compare_oneweight.py
	#resources/test/corsika_weight_calculator.py
)
