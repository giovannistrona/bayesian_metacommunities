# bayesian_metacommunities

Code accompanying the manuscript "A Bayesian network approach to trophic metacommunities shows habitat loss accelarates top species extinction rates", by Johanna Häussler, György Barabás and Anna Eklöf. 

The directory "c_code" contains the code in C++ and Bash needed to generate the food webs (binary adjacency matrices). To generate the webs, run the "build_webs.sh" file. Successfully running this file requires:

* Installation of C++
* The gsl library

The directory "code" contains all the source code in R and C++ needed to replicate the results in the manuscript. 
To run the simulations, run the "launch.R" file. This will create all necessary data and run all additional code automatically. Among others, it calls the function "run_instance.R" in which all the magic happens. Successfully running this file requires:

* An installation of R (version 3.0 or higher should work)
* Four additional R packages (igraph; NetIndices; Rcpp; tidyverse)

To generate the figures, run the "figure_script_manuscript.R" script. This file contains all code needed to recreate the figures in the manuscript; "figure_script_SI.R" contains code examples to recreate some figures in the SI. Successfully running these files requires:

* An installation of R (version 3.0 or higher should work)
* Four additional R packages (cowplot; magick; RColorBrewer; tidyverse)

The directory "code/ryser_code" contains the source code in R and C++ needed to replicate the results in the SI, S8. 
To run the simulations, run the "launch.R" file. This will create all necessary data and run all additional code automatically. Among others, it calls the script "run_instance.R" in which all the magic happens. Successfully running this file requires:

* An installation of R (version 3.0 or higher should work)
* Four additional R packages (igraph; NetIndices; Rcpp; tidyverse)

To generate the corresponding figures in the SI (Fig. S16, Fig. S17), use the file "plot_heatmaps.R". This requires the attached package scales to be at least of version 1.0.1 and the farver package to be 2.0.3.  

The directory "code/link_removal" contains the source code in R and C++ needed to replicate the results in the SI, S9. To run the simulations, run the "launch.R" file. This will create all necessary data and run all additional code automatically. Successfully running this file requires:

* An installation of R (version 3.0 or higher should work)
* Four additional R packages (igraph; NetIndices; Rcpp; tidyverse)

To generate the corresponding figure in the SI (Fig. S18), use the file "figure_script_linkrem.R". 

The file "code/launch_gridlike.R" contains all the source code in R needed to replicate the results in the SI, S10. To run the simulations, run the "launch_gridlike.R" file. This will create all necessary data and run all additional code automatically. Among others, it calls the function "run_instance.R" in which all the magic happens. Successfully running this file requires:

* An installation of R (version 3.0 or higher should work)
* Four additional R packages (igraph; NetIndices; Rcpp; tidyverse)

To generate the corresponding figure in the SI (Fig. S19), use the file "analyze_gridlike.R". 

All included code is covered under the GPL, version 3, available here: (http://www.gnu.org/licenses/gpl-3.0.en.html).

Any question about the code can be sent to Johanna Häussler (johanna.haeussler@idiv.de) and György Barabás (gyorgy.barabas@liu.se). 
