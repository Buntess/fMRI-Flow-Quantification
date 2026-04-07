# fMRI Inflow Effect Quantification

This repository contains MATLAB tools and scripts for simulating and quantifying inflow effects in functional Magnetic Resonance Imaging (fMRI). The project models CSF flow signal in different geometries (pipes and funnels) and evaluates the robustness of these models against noise and geometric assumptions.

## Project Structure

The repository is organized into core simulation functions, signal generation tools, and scripts for validation and investigation.

### Core Configuration & Data
* `generate_data.m`: Initializes, creates, and stores the fundamental variables required for the simulations.
* `get_protocol_settings.m`: Returns a `struct` containing the MRI modeling parameters and protocol settings used throughout the project.

### Signal Simulation
* `get_signal_pipe.m`: Simulates fMRI inflow signals for a standard pipe (cylindrical) geometry.
* `get_signal_funnel.m`: Simulates fMRI inflow signals for a funnel (tapered) geometry.
* `get_signal_profile.m`: Generates the specific inflow signal corresponding to a given flow profile.

### Experimental Investigations
* `investigate_N.m`: A script that investigates the impact of varying spatial resolutions on the simulated data.
* `investigate_profile.m`: A script designed to analyze how different fluid flow profiles affect the resulting fMRI inflow signal.

## Validations

The repository includes two main validation pipelines to test the efficacy and limits of the models.

### 1. Pipe Validation (Noise Robustness)
**Files:** `pipe_validation.m`, `pipe_plots.m`
These scripts represent the primary workflow for validating the pipe geometry model. The main objective of this validation is to test the model's resilience and determine how robust the pipe modeling is when subjected to noise.

### 2. Funnel Validation (Geometric Approximation)
**Files:** `funnel_validation.m`, `funnel_plots.m`
These scripts handle the validations for the funnel geometry. The core purpose of this validation is to investigate whether a complex funnel geometry can be accurately and reliably simulated using the simpler pipe modeling approach.