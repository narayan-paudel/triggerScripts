#!/bin/bash

ebin=7.0
rad=$(python -c "print(800+(int($ebin)-5)*300 + (2*(int($ebin)-5)//3)*300 + ((int($ebin)-5)//3)*300)")
echo rad $rad