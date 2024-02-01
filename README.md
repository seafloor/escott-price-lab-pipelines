# Neuro-genetic pathways
Standardised pathway annotations for use in genetics of brain disorders.

## Overview

This repository is dedicated to the development and maintenance of a comprehensive pipeline for downloading, processing, and standardizing genetic pathway data. Our goal is to provide researchers and enthusiasts with easy access to genetic pathways in a standardized format, making it simpler to integrate and compare data across different studies and genome builds.

The pipeline is designed to automate the liftover process to different genome builds, ensuring that the data is accessible and usable regardless of the reference genome version. Additionally, we provide versions of the data with and without the APOE region, catering to various research needs and preferences.

**Please note:** This project is currently a work in progress. We are actively developing and refining the pipeline and its documentation. Feedback and contributions are highly welcomed!

## Features

- **Automated Downloading:** Scripts to automatically fetch genetic pathway data from multiple sources.
- **Standardized Formatting:** Conversion of downloaded data into a standardized format for uniform access and analysis.
- **Genome Build Liftover:** Automated liftover processes to convert data between different genome builds.
- **APOE Region Handling:** Provision of data versions with and without the APOE region to support diverse research objectives.
- **Documentation:** Detailed documentation for each step of the pipeline, ensuring clarity and ease of use.

## Getting Started

To get started with the Genetic Pathways Standardization Pipeline, please follow the steps below:

1. **Clone the Repository:**

```bash
git clone https://github.com/your-repo/genetic-pathways-standardization.git
cd genetic-pathways-standardization
```

2. **Install Dependencies:**

Make sure you have Python 3.x installed and then install the required Python libraries:

```bash
pip install -r requirements.txt
```

3. **Configure Settings:**

Edit the `config.json` file to specify the sources from which to download the data and the target genome builds for the liftover process.

4. **Run the Pipeline:**

Execute the main script to start the downloading and processing of genetic pathways:

```bash
python run_pipeline.py
```

## Contributing

We welcome contributions 💚 Whether you have suggestions for new features, improvements, or bug fixes, your input is valuable to us.

## License

This project is licensed under the GPL v.3 licence - see the [LICENSE](LICENSE) file for details.
