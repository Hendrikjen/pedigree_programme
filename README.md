# pedigree programme

[xxx] is a pedigree analysis tool which was developed and implemented as part of a bioinformatic's master thesis in 2023 at Leipzig University. It is a C++ written console application which was designed to calculate dyadic relatedness coefficients from a given pedigree without being limited by the number of considered generations, the number of individuals, or the depth and incompleteness of the pedigree itself. Additionally, the programme provides some further informations about the respective relatedness paths between the focal individuals, for instance the name and sex of ancestors along the path, the lowest common ancestors, the kinlabel, or the minimal  detectable inbreeding value for each individual. The functionality and accuracy was adequately tested with multiple simulated populations as well as with a real pedigree of over 12 000 rhesus macaques from Cayo Santiago.

Since behavioural ecologists often have to deal with incomplete pedigree, due to unknown sires, the second part of the programme focus on an implementation of an adapted simulated annealing algorithm to find the best fully-reconstructed pedigree solution based on realized relatedness values obtained from whole genome sequencing. Eventually, it aims to provide a pedigree without gaps for which the difference between the provided realized relatedness values and the simultaneously calculated pedigree-derived relatedness coeffiecients is minimal over all dyads (see more informations in the section _Implementation/Simulated annealing_)

## Getting Started
<details>
<summary>Installation guide</summary>
  
- to get the programme, download (and don't forget to unzip) the repository to your local filesystem
- after downloading the source code, open the command line and navigate within the terminal into the folder _pedigree_programme/source/_
  - you can check with `ls` if you are in the correct folder if there are multiple Headerfiles (.h) and the respective source code files (.cpp) as well as _main.cpp_ and the makefile _makefile_pedigree_programme_
- run in the command line `make -f makefile_pedigree_programme`
- use the command `./pedigree_programme` to start the programme (depending on what you want to do, you have to add further arguments after the command)
- for general information you can type `./pedigree_programme -h`to list all possible command line arguments, or `./pedigree_programme -v` to get the current version
</details>

<details>
<summary>Command line arguments</summary>

`-f <functionality>` [string]
  - **options**: 
    - _relatedness_: calculates the dyadic relatedness (+ path characteristics) from a given pedigree
    - _simulation_: simulates a pedigree
    - _annealing_: starts a simulated annealing algorithm to fill the parental gaps within a pedigree based on realized relatedness values
  - **default**: [empty] (the programme starts without task)

<details><summary>
functionality == relatedness</summary>

#### required arguments
- `-p <input_pedigree>` [string]: path to pedigree file, e.g. _pedigree.txt_

#### optional arguments
- `-c <cores>` [int]
  - **options**: number of cores for multiprocessing
  - **default**: 1 (no multiprocessing)
- `-d <input_dyadlist>` [string]
  - **options**: path to file with selected dyads e.g. _dyad_selection.txt_
  - **default**: [empty] (all dyads within the pedigree will be analysed)
- `-e <output_extend>` [string]
  - **options**:
    - _full_: returns the full dyadlist output, including path characteristics
    - _reduced_: returns only dyadlist with dyadic relatedness coefficients
  - **default**: full
- `-l <generation_limit>` [int]
  - **options**: restricts the distance to potential lowest common ancestors, e.g. if generation_limit == _3_, only paths up to the grandparent generation will be returned, great-grand-parents will be considered as unrelated
  - **default**: [empty] (no limitation; all ancestors of a focal will be considered as potential lowest common ancestor)
- `-o <output>` [string]
  - **options**: custom output name (prefix) e.g. if output == _programme_output_, the resulting output files will be named "programme_output_dyadlist.txt" and "programme_output_info.txt"
  - **default**: [empty] (the input file name will be used as prefix)
- `-r <reduce_node_space>` [bool]
  - **options**: 
    - _true_: before calculating the dyadic relatedness, the number of individuals will be reduced which means that only descendants of the focal's common ancestors will be considered in the analysis (it effectively reduces the search space without affecting the result, but might be only beneficial in almost completely reconstructed pedigrees with a long history due to the extra computational cost)
    - _false_: no prior narrowing of the search space
  - **default**: false

#### Example
`./pedigree_programme -f relatedness -p pedigree.txt -e reduced -c 5`

</details>
<details><summary>
functionality == simulation</summary>

#### required arguments
- `-n <start_individual>` [int]: number of individuals at the start of the simulation
- `-s <simulation_duration>` [int]: number of years to restrict the duration of the simulation

#### optional arguments
- `-a <max_age>` [int]
  - **options**: age maximum in population (individuals who reach the maximum age will decease in the following year)
  - **default**: 200
- `-b <birth_rate>` [double]
  - **options**: specifies the annual increment in the number of offsprings born each year during the population simulation
  - **default**: 4.0
- `-q <death_rate>` [double]
  - **options**: specifies the annual increment in the number of deaths each year during the population simulation
  - **default**: 3.0
- `-y <default_year>` [int]
  - **options**: start year for population simulation
  - **default**: 1900

#### Example
`./pedigree_programme -f simulation -n 20 -s 10 -y 1938`
</details>
<details><summary>
functionality == annealing</summary>

#### required arguments
- `-d <input_dyads_complete>` [string]: path to dyadlist with realized relatedness values, e.g. _true_dyads.txt_
- `-p <input_pedigree>` [string]: path to pedigree file (with gaps), e.g. _pedigree.txt_

#### optional arguments
- `-c <cores>` [int]
  - **options**: number of cores for multiprocessing
  - **default**: 1 (no multiprocessing)
- `-i <init_temp>` [double]
  - **options**: start temperature 
  - **default**: [empty] (automatically calculated)
- `-t <stop_temp>` [double]
  - **options**: stop temperature, if current temperature falls below stop temperature, the algorithm terminates
  - **default**: 1.0
- `-x <temp_decay>` [double]
  - **options**: the temperature multiplication factor to determine the number of iterations (if the number of iteration _n_ is set, the decay factor can be calculated with temp_decay = $\sqrt[n]{\frac{t_{stop}}{t_{init}}} $
  - **default**: 0.99

#### Example
`./pedigree_programme -f annealing -p pedigree_with_gaps.txt -d realized_dyadic_relatedness.txt -x 0.995 -c 5 -m 1000 -w 1000`
</details>

#### general optional arguments
- `-g <gestation_length>` [int]
  - **options**: gestation length in days
  - **default**: 200
- `-m <maturation_age_m>` [int]
  - **options**: maturation age of males in days
  - **default**: 1250
- `-w <maturation_age_f>` [int]
  - **options**: maturation age of females in days
  - **default**: 1095

</details>


## Input requirements
<details>
<summary>Relatedness calculation</summary>

#### Pedigree
- Input file format: .txt (tab-separated)
- no header
- empty NA values (like "") lead to adverse behaviour or programme abort
- columns (order and format is mandatory): ID, sex, birthseason/year, mom_ID, sire_ID, day of birth (DOB), day of death (DOD), nonsire, nondam

|column|data type|missing value|comment|
|-|-|-|-|
|ID |string| ID names like _UNK_, _NA_, _unknown_, _unkn_f_, and _unkn_m_ have to be avoided|ID names have to be unique and have to be unambiguously assignable to parent IDs; every parent ID from mom_ID or sire_ID has to be listed in the pedigree separately
|sex |char| u| usage of the following options only _f_ = female, _m_ = male, or _u_ = unknown sex
|birthseason |int|0| year
|mom_ID |string|unknown| have to be relatable to exactly one ID, respectively one female individual in the pedigree
|sire_ID |string|unknown| have to be relatable to exactly one ID, respectively one male individual in the pedigree
|DOB |string (dateformat)| NA| in the format: 01-01-1900
|DOD |string (dateformat)|NA| in the format: 01-01-1900
|nonsire |string| NA|IDs of previously excluded sires strung together (have to be relatable to exactly one ID of the respective sex in the pedigree); separated by @ e.g. _indiv1@indiv2@indiv3_; ensure that each individual has at least one remaining potential sire within the pedigree, else the individual will be assumed to be a founder individual
|nondam |string| NA|IDs of previously excluded moms strung together (have to be relatable to exactly one ID of the respective sex in the pedigree); separated by @ e.g. _indiv1@indiv2@indiv3_; ensure that each individual has at least one remaining potential mother within the pedigree, else the individual will be assumed to be a founder individual|

- [example](example/example_input_pedigree.txt)

#### Dyad Selection
- Input file format: .txt (tab-separated) 
- no header
- empty NA values (like "") lead to adverse behaviour or programme abort
- columns (order and format is mandatory): ID_1, ID_2
  - ID names have to be unique and have to be unambiguously assignable to pedigree IDs; every focal ID has to be listed in the pedigree separately; ID names like _UNK_, _NA_, _unknown_, _unkn_f_, and _unkn_m_ have to be avoided
- [example](example/example_input_dyad_selection.txt)
</details>
<details>
<summary>Simulated Annealing</summary>
  ...
</details>

## Example
<details>
<summary>Relatedness calculation</summary>
<p align="center">
  <img src="example/mini_example_git.png" width="300">
</p>
<details>
<summary>Input/Output files</summary>

<details>
<summary> Input file (pedigree)
</summary>

|ID|sex|birthseason|mom|sire|DOB|DOD|nonsire|nondam|
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|A|f|1905|unknown|unknown|01-01-1900|NA|NA|NA|
|B|f|1911|A|unknown|01-01-1911|NA|NA|NA|
|C|m|1912|unknown|unknown|01-01-1912|NA|NA|NA|
|D|f|1913|A|unknown|01-01-1913|NA|NA|NA|
|E|f|1914|A|unknown|01-01-1914|NA|NA|NA|
|F|m|1915|unknown|unknown|01-01-1915|NA|NA|NA|
|G|m|1920|B|unknown|01-01-1920|NA|NA|NA|
|H|f|1921|D|C|01-01-1921|NA|NA|NA|
|I|m|1922|E|F|01-01-1922|NA|NA|NA|
|J|m|1923|E|F|01-01-1923|NA|NA|NA|
|K|m|1928|H|G|01-01-1928|NA|NA|NA|
|L|f|1929|H|I|01-01-1929|NA|NA|NA|

[example_input_pedigree.txt](example/example_input_pedigree.txt) 
</details>

<details>
<summary> Input file (dyad selection)
</summary>

|ID_1|ID_2|
| ------------- | ------------- |
|C|F|
|H|L|
|I|J|
|K|L|
|C|G|
|D|G|
|D|J|

[example_input_dyad_selection.txt](example/example_input_dyad_selection.txt) 
</details>

<details>
<summary> Output file (pedigree): additional pedigree info like generational depth and minimal inbreeding value
</summary>

|ID|sex|BS|mom|sire|DOB|DOD|pot_sire|pot_mom|full_generations|min_f|
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|A|f|1905|unkn_f|unkn_m|1-1-1900|0-0-0|NA|NA|1|0.000000000000000|
|B|f|1911|A|unkn_m|1-1-1911|0-0-0|NA|NA|1|0.000000000000000|
|C|m|1912|unkn_f|unkn_m|1-1-1912|0-0-0|NA|NA|1|0.000000000000000|
|D|f|1913|A|unkn_m|1-1-1913|0-0-0|NA|NA|1|0.000000000000000|
|E|f|1914|A|unkn_m|1-1-1914|0-0-0|NA|NA|1|0.000000000000000|
|F|m|1915|unkn_f|unkn_m|1-1-1915|0-0-0|NA|NA|1|0.000000000000000|
|G|m|1920|B|unkn_m|1-1-1920|0-0-0|NA|NA|1|0.000000000000000|
|H|f|1921|D|C|1-1-1921|0-0-0|NA|NA|2|0.000000000000000|
|I|m|1922|E|F|1-1-1922|0-0-0|NA|NA|2|0.000000000000000|
|J|m|1923|E|F|1-1-1923|0-0-0|NA|NA|2|0.000000000000000|
|K|m|1928|H|G|1-1-1928|0-0-0|NA|NA|2|0.031250000000000|
|L|f|1929|H|I|1-1-1929|0-0-0|NA|NA|3|0.031250000000000|

[example_output_pedigree_info.txt](example/example_output_pedigree_info.txt)
</details>
<details>
<summary> Output file (dyadlist): path characteristics
</summary>

|ID 1|ID 2|dyad|relatedness coefficient|paths|pathline|kinline|LCA|depth|kinlabel|fullhalf|min_DGD|
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
|C|F|C_F|0|NA|NA|NA|NA|NA|nonkin|NA|1|
|H|L|H_L|0.531250000000000|H@L/@/H@D@A@E@I@L|ff/@/ffffmf|mat/@/mixed|H/@/A|0/1/@/2/3|daughter&mother/@/1st-cousins-once-removed|half/@/half|2|
|I|J|I_J|0.500000000000000|I@E@J/@/I@F@J|mfm/@/mmm|mat/@/pat|E/@/F|1/1/@/1/1|brothers/@/brothers|full/@/full|2|
|K|L|K_L|0.296875000000000|K@H@L/@/K@H@D@A@E@I@L/@/K@G@B@A@D@H@L/@/K@G@B@A@E@I@L|mff/@/mffffmf/@/mmfffff/@/mmfffmf|mat/@/mixed/@/mixed/@/mixed|H/@/A/@/A/@/A|1/1/@/3/3/@/3/3/@/3/3|siblings/@/2nd-cousins/@/2nd-cousins/@/2nd-cousins|half/@/half/@/half/@/half|2|
|C|G|C_G|0|NA|NA|NA|NA|NA|nonkin|NA|1|
|D|G|D_G|0.125000000000000|D@A@B@G|fffm|mat|A|1/2|nephew&aunt|half|1|
|D|J|D_J|0.125000000000000|D@A@E@J|fffm|mat|A|1/2|nephew&aunt|half|1|

[example_output_dyadlist.txt](example/example_output_dyadlist.txt)
</details>
</details>

<details>
<summary>Output explanation (path characteristics)</summary>

To further explain the column in the dyadlist output, we will look on the examplary dyad (E_G) from the pedigree example above. The focal individuals E (circle = female) and G (square = male) are related only by maternal ancestors (kinline = mat), whereby the lowest common ancestor A is one edge apart from E and two from G (depth = 1/2) which codes in combination with the sex for the kinlabel nephew/aunt. Each focal has at least one unknown parent, therefore the min DGD is 1.

|name | explanation | example |
| ------------- | ------------- | ------------- |
|path | consecutive list of nodes along the relatedness path (edge directions are left unregarded) | E@A@B@G|
|lca | lowest common ancestor within path | A |
|pathline | sequence of sexes (f/m/u) along the path | fffm |
|kinline | whether the path consists solely of maternal or paternal ancestors; “mixed” if both sexes occur | mat |
|depth | path length from LCA to each focal | 1/2 |
|kinlabel | kinclass label based on the table of consanguinity (see below) | nephew-aunt |
|fullhalf | whether two identical paths exist with different common ancestors, e.g. differentiation between full- and half-siblings | half |
|min\_DGD | minimal dyadic genealogical depth states the pedigree completeness for the dyad; i.e. the minimal amount of fully resolved generations starting from both focals | 1 |


<p align="center">
<img src="https://upload.wikimedia.org/wikipedia/commons/0/0d/Table_of_Consanguinity_showing_degrees_of_relationship.svg" width="500">

https://upload.wikimedia.org/wikipedia/commons/0/0d/Table_of_Consanguinity_showing_degrees_of_relationship.svg
</p>
  
</details>
</details>
<details>
<summary>Population Simulation</summary>

examplary output of a simulated pedigree with 20 founder individuals born/start in 1950, simulated for 10 years: [simulated pedigree](example/population_simulation/example_simulation.txt) and the respective list of [dyadic relatedness coefficients](example/population_simulation/example_simulation_dyadic_paths.txt) for all 1442 individuals.
</details>
<details>
<summary>Simulated Annealing</summary>
  ...
</details>

## Implementation
<details>
<summary>Relatedness Coefficient</summary>

#### Recursive relatedness coefficient calculation [^1]

To calculate the dyadic relatedness coefficient, the pedigree G is conceived as as a directed, acyclic graph, consisting of two distinct classes of vertices, $V_1$ (males) and $V_2$ (females) whereas each vertex represents an individual. Edges within the graph referred to one-directional, direct kinship bonds between parent and offspring, which implies that for each (heterogamous) node at least two edges exist (to the mother and to the father), or more in case of own offspring. But while in reality, pedigrees often consists of missing parents, two imaginary nodes $\rho_1\ \epsilon\ V_1$ and $\rho_2\ \epsilon\ V_2$ are added, serving as a compensatory substitute for unknown mothers or sires.

Generally, the relatedness coefficient of an individual $x\ \epsilon\ V$ to itself is stated as $f\left(x,x\right)=1$ while the relatedness of two different focals $f\left(x,y\right)$ can be expressed by the following recursive formula 
$$f\left(x,y\right)=\ \frac{1}{4}\left[f\left(x_1,y_1\right)+f\left(x_1,y_2\right)+f\left(x_2,y_1\right)+f(x_2,y_2)\right]$$ ($x_1,\ x_2$ as parents of $x$; $y_1,\ y_2$ as parents of $y$ while $x_1,\ y_1\ \epsilon\ V_1$ and $x_2,\ y_2\ \epsilon\ V_2$). 
In the particular case of determining the relatedness coefficient between an individual $x$ and its ancestor $x_i$, it is calculated by
$$f\left(x,x_i\right)=\ \frac{1}{2}\left[f\left(x_1,x_i\right)+f\left(x_2,x_i\right)\right]$$
($x,\ x_i\ \epsilon\ V;\ x_1\ \epsilon\ V_1$ and $x_2\ \epsilon\ V_2$ as parents of $x$). Even more specific, if $x_i \equiv x_1 \lor x_2$, the relatedness between parent and offspring is given by 
$$f\left(x,x_1\right)=\ \frac{1}{2}\left[1+f\left(x_1,x_2\right)\right]$$
At last, in case of imaginary nodes, $\rho_1$ and $\rho_2$ are assumed as unrelated to each other or any other individual $x\ \epsilon\ V:$ 
$$f\left(\rho_1,\rho_2\right)=f\left(x,\rho_1\right)=f\left(x,\rho_2\right)=0$$
Based on these functions, the programme computes the relatedness between a dyad step by step until it either identifies their lowest common ancestor or terminates due to a trivial solution.
</details>

<details>
<summary>Simulated Annealing </summary>

#### Adapted Simulated Annealing Algorithm [^1]

Within the programme a simulated annealing algorithm is implemented to fill possibly existing gaps within a given pedigree. Therefor, it uses the discrepancy between user-provided realized relatedness values (for instance obtained from whole genome sequencing) and the calculated pedigree-derived relatedness values as cost function. While trying to minimize the cost/discrepancy by simulated annealing, the aim is to find the pedigree solution which explains best the variance, especially in case of missing ancestors which can be accompanied with an underestimation of relatedness values. 
$$F =\Sigma\ |f(x,y) - g(x,y) | \to min$$ (with $f(x,y)$ as the pedigree-based dyadic relatedness and $g(x,y)$ as the dyadic realized relatedness)
To fit the specific problem, the general simulated annealing algorithm is adapted as explained in the following outline:

- At first, all pedigree gaps need to be identified.
- Create an initial solution by randomly assigning parents from a pool of suitable candidates for each gap. Suitable candidate are parents who were alive and mature at the time of conception (male) or birth (female) and were not excluded as potential parent priorly due to genetic analysis or because a female has already an offspring in the respective cohort.
- Calculation of the relatedness coefficient for each relevant dyad (those for which realized relatedness values are available)
- Evaluate the difference between the realized and pedigree-derived relatedness values of the initial solution for each relevant dyad
- Save the current difference as the best-known difference, and the initial solution as the best pedigree.
- Iteration: While the current temperature is above the (given) stop temperature:
  - Create a new solution by exchanging one potential parent with another suitable candidate.
  - Calculate new relatedness values for dyads affected by this change (all relevant dyads which include the offspring, the old and the new parent candidate).
  - Compare the old and new relatedness values to determine the discrepancy between the current and the new solution.
  - If the new solution is worse, apply the Metropolis acceptance criterion to decide whether to accept it or not: $$e^\frac{F_{n}-F_{c}}{T} > X\to [0,1]$$ 
 (with $F_n$ as fitness function of the new solution and $F_c$ of the current solution; $T$ as temperature and $X$ as random number in the range between 0 and 1)
  - If accepted (or the new solution is better in the first place), the new solution becomes the current solution; otherwise, it's rejected, and the previous solution remains in place.
  - If necessary, update the best difference and pedigree.
- Finally, save the last pedigree solution in a file.


</details>

## Contribution and citation
I want to thank and express my deep gratitude to all people who contributed to this project and who never failed to come up with helpful recommandations, fresh ideas and new perspectives. Especially, the unconditional support, encouragement and guidance of my supervisors A. Freudiger, T. Gatter, P. Stadler and A. Widdig was essential for bringing this project to a successful end. 

Please use the BibTex format, provided by Github or cite this programme as 

**Westphal, H. (2023). Pedigree programme (Version 1.0.0) [Computer software].** _https://github.com/Hendrikjen/pedigree_programme_

Further background information may be found in my [master thesis](master%20thesis/MA_Assessing-dyadic-relatedness-in-rhesus-macaques-using-pedigree-data_HWestphal_2023.pdf).

Contact email: hw53vake@studserv.uni-leipzig.de

[^1]: Westphal, H. (2023). Assessing dyadic relatedness in rhesus macaques using pedigree data [Master thesis]. [_https://github.com/Hendrikjen/pedigree_programme/master thesis/MA_Assessing-dyadic-relatedness-in-rhesus-macaques-using-pedigree-data_HWestphal_2023.pdf_](master%20thesis/MA_Assessing-dyadic-relatedness-in-rhesus-macaques-using-pedigree-data_HWestphal_2023.pdf)
 
