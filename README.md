# Pedigree Analysis Framework

A Python-based framework for analyzing pedigrees and estimating inheritance patterns using minimal phenotypic data.

## Features

- Analyze pedigrees using only health status information (affected/unaffected)
- Estimate inheritance patterns (Autosomal Dominant, Autosomal Recessive, Y-linked)
- Calculate p-values using Monte Carlo simulation
- Parallel processing for improved performance on large pedigrees
- Simple .ped file input format

## Installation

```bash
git clone https://github.com/yourusername/pedigree-analysis-framework.git
cd pedigree-analysis-framework
pip install -r requirements.txt
```

## Usage

```python
from pedigree import Pedigree, Models

# Load pedigree data
ped = Pedigree()
ped.read_from_file("path/to/your/pedigree.ped")

# Analyze for Autosomal Dominant inheritance
p_value = ped.calc_p_value(Models.AutosomalDominant)
print(f"P-value for Autosomal Dominant inheritance: {p_value}")

# Analyze for Autosomal Recessive inheritance
p_value = ped.calc_p_value(Models.AutosomalRecessive)
print(f"P-value for Autosomal Recessive inheritance: {p_value}")
```

## Input File Format

The framework accepts .ped files with the following format:

```
ID,Gender,Status,Father,Mother
0,M,U,N,N
1,F,A,N,N
2,M,U,0,1
3,F,A,0,1
```

- ID: Unique identifier for each individual
- Gender: M (Male) or F (Female)
- Status: A (Affected) or U (Unaffected)
- Father: ID of the father (N if unknown)
- Mother: ID of the mother (N if unknown)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use this framework in your research, please cite:

[Your Name]. (2024). Pedigree Analysis Framework. GitHub repository: https://github.com/yourusername/pedigree-analysis-framework

## Contact

[Your Name] - [your.email@example.com]

Project Link: [https://github.com/yourusername/pedigree-analysis-framework](https://github.com/yourusername/pedigree-analysis-framework)
