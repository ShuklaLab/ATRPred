# ATRPred: <ins>A</ins>nti-TNF <ins>T</ins>reatment <ins>R</ins>esponse <ins>Pred</ins>ictor

A package for predicting whether a Rheumatoid Arthritis patient will respond to anti-TNF treatment.

**Visit us at:** [Shukla Lab](https://shuklalab.github.io/)

*Software author/maintainer: [@bodhayan](https://github.com/bodhayan)*

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

## Contributing and licensing
[MIT](https://choosealicense.com/licenses/mit/) License is free license software and grants the software end user rights such as copying, modifying, merging, distributing, etc.

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Publication
If you're using our tool and/or dataset, please cite:

[Prasad B, McGeough C, Eakin A, Ahmed T, Small D, Gardiner P, Pendleton A, Wright G, Bjourson AJ, Gibson DS, Shukla P. ATRPred: A machine learning based tool for clinical decision making of anti-TNF treatment in rheumatoid arthritis patients. PLoS Comput Biol. 2022 Jul 5;18(7):e1010204. doi: 10.1371/journal.pcbi.1010204. PMID: 35788746; PMCID: PMC9321399.](https://doi.org/10.1371/journal.pcbi.1010204)
***
*Last updated on: 6 Sep 2022*
