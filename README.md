# Bayesian predictive model for the lifetime of indoor ultrafine particle

This repository contains MATLAB code tailored for Bayesian parameter identification in indoor ultrafine particle modeling. Leveraging advanced Bayesian inference and detailed particle dynamics models, this toolkit enables robust analyses and predictions of indoor ultrafine particle behaviors.

## Detailed Description

### 1. Initial Settings (`BPM_IndoorUFP_InitialSetting`)
The `BPM_IndoorUFP_InitialSetting` function sets initial parameters for ultrafine particle modeling based on particle sources, integrating experimental data to establish baseline configurations for simulation.

### 2. Data Handling and Bayesian Setup (`BPM_IndoorUFP_Main.m`)
The `BPM_IndoorUFP_Main.m` module manages experimental data, setting up a structured Bayesian analysis framework tailored for ultrafine particle dynamics.

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
git clone https://github.com/MPMC-Lab/BPM_IndoorUFP.git
```
Alternatively, the source files can be downloaded through github menu 'Download ZIP'.

## Contributing

Contributions are welcome. Please fork the repository and submit a pull request for review.

## License

This project is under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## References

For further information, refer to associated research papers and visit the School of Mathematics and Computing at Yonsei University.

## Contact Information

Should you require assistance or have any inquiries, kindly contact Dr. Jung-Il Choi via email at [jic@yonsei.ac.kr](mailto:jic@yonsei.ac.kr).

For detailed information and further reading related to BPM_GasAdsorption, you are encouraged to consult the reference paper. Additional insights and resources are available at School of Mathematics and Computing, within the domain of Computational Science and Engineering at Yonsei University. For more details, please visit our website: [mpmc.yonsei.ac.kr](http://mpmc.yonsei.ac.kr).

