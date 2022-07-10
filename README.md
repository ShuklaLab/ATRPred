# ATRPred: <ins>A</ins>nti-TNF <ins>T</ins>reatment <ins>R</ins>esponse <ins>Pred</ins>ictor

A package for predicting whether a Rheumatoid Arthritis patient will respond to anti-TNF treatment.

**Visit us at:** [Shukla Lab](https://shuklalab.github.io/)

![GitHub](https://img.shields.io/github/license/ShuklaLab/ATRPred)
![GitHub R package version](https://img.shields.io/github/r-package/v/ShuklaLab/ATRPred)

## Requirements

ATRPred was coded in R v3.6, and require a couple of packages (`caret` and `RANN`) to be pre-installed, which can be installed as below:

```
install.packages(c("caret","RANN"))
```

## Installation and Updation

To install the latest version of ATRPred, please follow the steps:

```
install.packages("devtools")
devtools::install_github("ShuklaLab/ATRPred")
```

To update the package just run the same command again.

## Running the function

To run the ATRPred, please run the following code (input template along with sample inputs are present in examples folder):
```
library("ATRPred")
antiTNFresponse()
```

Or directly, by typing:

```
ATRPred::antiTNFresponse()
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[MIT](https://choosealicense.com/licenses/mit/)


## Publication
[Prasad B, McGeough C, Eakin A, Ahmed T, Small D, Gardiner P, Pendleton A, Wright G, Bjourson AJ, Gibson DS, Shukla P. ATRPred: A machine learning based tool for clinical decision making of anti-TNF treatment in rheumatoid arthritis patients. PLoS Comput Biol. 2022 Jul 5;18(7):e1010204. doi: 10.1371/journal.pcbi.1010204. Epub ahead of print. PMID: 35788746.](https://doi.org/10.1371/journal.pcbi.1010204)
