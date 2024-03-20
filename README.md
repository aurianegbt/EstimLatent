# Enhancing mechanistic modeling using longitudinal transcriptomics data

We generated $R = 100$ datasets mimicking the trajectories of $S$ cells and antibodies using our mechanistic model in 50 or 10 individuals (see simulationGroup1 and simulationGroup2 in simulx project). Sampling time were taken at  at days $\{0, 7, 14, 21, 28, 35\}$ for antibodies and at days $0$ to $14$ and $21$ to $35$ for $S$ cells. The true parameter values $\theta^*$ for the mechanistic model are detailed below. We explored three distinct scenarios: one where information is solely derived from antibody data, and two where information is derived from both antibody data and transcriptomic markers after deconvolution, with varying levels of noise ($G_1$ for the high noise setting, $G_2$ for low noise setting, either made by choosing different $\alpha$ or different levels of noise $\sigma^2$. 

# Simulation settings 

## For simulation named "Mult" : 

```r
print(pop = read.csv("Model/IreneModelMult/Simulation/populationParameters.txt")[1,2:10,drop=F])
```

## For simulation named "MP" : 

```r
print(pop = read.csv("Model/IreneModelMP/Simulation/populationParameters.txt")[1,2:10,drop=F])
```

## For simulation named "Comp" : 

```r
print(pop = read.csv("Model/IreneModelComp/Simulation/populationParameters.txt")[1,2:10,drop=F])
```
