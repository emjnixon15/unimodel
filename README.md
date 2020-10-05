# Stochastic transmission model for COVID-19 transmission in a university settings

## About the project
This code accompanies the paper "High COVID-19 transmission potential associated with re-opening universities can be mitigated with layered interventions" by Ellen Brooks-Pollock, Hannah Christensen, Adam Trickey, Gibran Hemani, Emily Nixon, Amy Thomas, Katy Turner,  Adam Finn, Matt Hickman, Caroline Relton, Leon Danon, available at https://doi.org/10.1101/2020.09.10.20189696. 


## Running the code 
The code is written in R version 4.0.2 and C. Steps for 

* Clone the repository 
```ShellSession
[user@host ~]$ git clone https://github.com/ellen-is/unimodel.git
```
* Open "UniModel.R" in R or Rstudio. You will need to install the R libraries Rcpp, tidyverse and gplots. 
* Run script. 

## Demo run and output
The code runs ten iterations of the model with default parameters and plots the output in around 2 seconds. Increase the number of simulations with the parameter nsim, found in UniModel.R. 
```ShellSession
nsim=100
```

The left panel shows the number of asymptomatic and symptomatic cases over time, with one line per trajectory. The right panel shows the average (of nsim runs) number of symptomatic cases per year group. 

## License
GPL3

## Contact
Ellen.Brooks-Pollock@bristol.ac.uk
