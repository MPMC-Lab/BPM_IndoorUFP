# Bayesian predictive model for for Indoor Ultrafine Particle Modeling

This repository contains MATLAB code tailored for Bayesian parameter identification in indoor ultrafine particle modeling. Leveraging advanced Bayesian inference and detailed particle dynamics models, this toolkit enables robust analyses and predictions of indoor ultrafine particle behaviors.

## Detailed Description

### 1. Initial Settings (`PI_IndoorUFP_InitialSetting`)
The `PI_IndoorUFP_InitialSetting` function sets initial parameters for ultrafine particle modeling based on particle sources, integrating experimental data to establish baseline configurations for simulation.

### 2. Data Handling and Bayesian Setup (`PI_IndoorUFP_Main.m`)
The `PI_IndoorUFP_Main.m` module manages experimental data, setting up a structured Bayesian analysis framework tailored for ultrafine particle dynamics.

### 3. IUQ Method Implementation (`uq_UFP`)
This function employs the Inverse Uncertainty Quantification method for accurate parameter estimation, processing multiple datasets for comprehensive analysis.

### 4. Forward Modeling (`Aerosol_Model`)
The `Aerosol_Model` function simulates ultrafine particle dynamics, using various models and solvers to predict particle behavior in indoor environments.

### 5. Supporting Functions
Essential computational functions support the modeling process, including particle distribution analysis and environmental impact calculations.

## Authors

- **Yesol Hyun** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - yesol2@yonsei.ac.kr
- **Jung-Il Choi** - School of Mathematics and Computing (Computational Science and Engineering), Yonsei University - jic@yonsei.ac.kr

## Installation

```bash
git clone https://github.com/MPMC-Lab/IndoorUFP_PI_repo.git
```
Alternatively, the source files can be downloaded through github menu 'Download ZIP'.

## Contributing

Contributions are welcome. Please fork the repository and submit a pull request for review.

## License

This project is under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

Thanks to all contributors and researchers for their valuable insights and feedback.

## References

For further information, refer to associated research papers and visit the School of Mathematics and Computing at Yonsei University.
