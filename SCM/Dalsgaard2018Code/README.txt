To run the extinction models, first download the Stochastic Coextinction Model code from the Supplementary Information of Vieira and Almeida-Neto (2015) (see references). The script netcascade (April 2014).R is required for the analyses in this paper.

The supplementary information contains three scripts:
1. Main.R – This script runs simulations on the 8 networks to calculate PE, SE and VE for each species in each network.
2. IterativeNodeDeletion.R – A wrapper around netcascade (April 2014).R that calls netcascade() repeatedly to iteratively remove species from the network.
3. IterNodeDelMultiSim.R – A wrapper around IterativeNodeDeletion.R that runs IterativeNodeDeletion() multiple times to account for the stochasticity in the netcascade() algorithm.

And four datasets:
1. data/rvalues.csv – The expert-assigned dependency categories for each plant species in each network used to determine R values
2. data/webs – Folder containing the 8 interaction matrices 
3. data/phylogeny - contains (1) a set of 1000 phylogenies for the 13 Caribbean hummingbird species that were analysed (Caribbean_Phylogenies_.nex), (2) the maximum clade credibility (MCC) phylogeny (Caribbean_MCC_Phylogeny.nex)
4. data/traits - species-level data for body mass and bill length

If you use the extinction modelling code please cite the following two papers: 
1. Dalsgaard et al. (in press) Trait evolution, resource specialisation and vulnerability to plant extinctions among Antillean hummingbirds. Proceedings of the Royal Society B.
2. Vieira, M.C. and Almeida‐Neto, M., 2015. A simple stochastic model for complex coextinctions in mutualistic networks: robustness decreases with connectance. Ecology Letters, 18(2), pp.144-152. doi: 10.1111/ele.12394

If you use the data please cite:
Dalsgaard et al. (in press) Trait evolution, resource specialisation and vulnerability to plant extinctions among Antillean hummingbirds. Proceedings of the Royal Society B.

References:
Vieira, M.C. and Almeida‐Neto, M., 2015. A simple stochastic model for complex coextinctions in mutualistic networks: robustness decreases with connectance. Ecology Letters, 18(2), pp.144-152. doi: 10.1111/ele.12394