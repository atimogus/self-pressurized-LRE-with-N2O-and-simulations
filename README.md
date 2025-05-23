## Overview

This repository contains code for the design and simulation of a self-pressurized liquid rocket engine. It automates the development of the entire propellant feed system, incorporating a concentric pressure vessel and a conical rocket nozzle.

The code supports both **cold flow** and **hot fire** simulations, enabling detailed performance evaluation and iteration during the design process.

## Usage

1. Open the `config.py` file and adjust the configuration variables according to your specific engine parameters.
2. At the bottom of the `config.py` script, you can specify the type of simulation you wish to perform.
3. To run the simulation, simply execute:

python main.py


## Features

- Automated propellant system design
- Support for concentric pressure vessel configurations
- Cold flow and hot fire simulation modes

## Roadmap / TODO

- [ ] Rocket trajectory analysis using **RocketPy**
- [ ] Detailed analysis of **showerhead (scrintle) injectors**
- [ ] Automated cooling channel design using **FreeCAD**
- [ ] Automated **CFD simulation** with **OpenFOAM**, including optimization of the number of cooling channels
