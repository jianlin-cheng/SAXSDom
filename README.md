# SAXSDom (updated by 2020/01/08)
This is a bioinformatics tool to use SAXS data to assemble protein domain structures into full-length structural models


**Installation**

**(1) Setup libraries (required)**

```
Updating

perl setup_env.pl

if report error as "Can't locate Env.pm in @INC (you may need to install the Env module)"

sudo dnf  install perl-Env
```
**(2) Configure MULTICOM system (required)**
```
perl configure.pl
```


**(3) Compile SAXSDom from source code (required)**
```
sh compile_SAXSDom.sh

#### output

- Boost version: 1.55.0
-- Boost version: 1.55.0
-- Found the following Boost libraries:
--   serialization
--   program_options
--   thread
-- CMAKE_CXX_FLAGS:         
running /data/jh7x3/SAXSDom/db_tools/tools/cmake-2.8.12.2/bin/cmake -E create_symlink "/data/jh7x3/SAXSDom/installation/Mocapy++-1.07/examples/data" "/data/jh7x3/SAXSDom/installation/Mocapy++-1.07/examples/data"  2>&1
-- Configuring done
-- Generating done
-- Build files have been written to: /data/jh7x3/SAXSDom/installation/Mocapy++-1.07
Scanning dependencies of target Mocapy
[  1%] Building CXX object src/CMakeFiles/Mocapy.dir/utils/LogFactorial.cpp.o
[  2%] Building CXX object src/CMakeFiles/Mocapy.dir/multinomial/MultinomialDensities.cpp.o
[  3%] Building CXX object src/CMakeFiles/Mocapy.dir/multinomial/MultinomialDensity.cpp.o
[  4%] Building CXX object src/CMakeFiles/Mocapy.dir/multinomial/MultinomialESS.cpp.o
[  6%] Building CXX object src/CMakeFiles/Mocapy.dir/dirichlet/DirichletDensities.cpp.o
[  7%] Building CXX object src/CMakeFiles/Mocapy.dir/dirichlet/DirichletDensity.cpp.o
[  8%] Building CXX object src/CMakeFiles/Mocapy.dir/dirichlet/DirichletESS.cpp.o
[  9%] Building CXX object src/CMakeFiles/Mocapy.dir/inference/abstractinfengine.cpp.o
[ 10%] Building CXX object src/CMakeFiles/Mocapy.dir/framework/dbn.cpp.o
...
...
...
Linking CXX static library ../libs/libMocapy.a
[ 96%] Built target Mocapy
Scanning dependencies of target SAXSDom
[ 97%] Building CXX object examples/CMakeFiles/SAXSDom.dir/SAXSDom.cpp.o
Linking CXX executable SAXSDom
[ 97%] Built target SAXSDom
Scanning dependencies of target mdarray
[ 98%] Building CXX object examples/CMakeFiles/mdarray.dir/mdarray.cpp.o
Linking CXX executable mdarray
[ 98%] Built target mdarray
Scanning dependencies of target mlr-uni
[100%] Building CXX object examples/CMakeFiles/mlr-uni.dir/mlr-uni.cpp.o
Linking CXX executable mlr-uni
[100%] Built target mlr-uni
```

**Example I**
```
cd examples

sh T1-run-3p02A.sh

#### output:

  All running jobs are done

  Ranking score for models can be found at /data/jh7x3/SAXSDom/test_out/3p02A_test//all_models_qprob/3p02A.Qprob_score

  Best model: /data/jh7x3/SAXSDom/test_out/3p02A_test//3p02A_SAXSDom_top1.pdb

Structure1: ./3p02A_SA  Length=  305
Structure2: ../../exam  Length=  305 (by which all scores are normalized)
Number of residues in common=  305
RMSD of  the common residues=    5.987

TM-score    = 0.6475  (d0= 6.41, TM10= 0.6060)
MaxSub-score= 0.4714  (d0= 3.50)
GDT-TS-score= 0.5795 %(d<1)=0.4557 %(d<2)=0.4754 %(d<4)=0.5377 %(d<8)=0.8492
GDT-HA-score= 0.4738 %(d<0.5)=0.4262 %(d<1)=0.4557 %(d<2)=0.4754 %(d<4)=0.5377
```



<h4> Run SAXSDom </h4>

```
Usage:
sh bin/run_SAXSDom.sh <target id> <fasta file> <saxs file> <domain list with path> <output directory> <number of cores to run, default:1>

Example:
sh bin/run_SAXSDom.sh 3p02A  examples/3p02A/3p02A.fasta  examples/3p02A/saxs_profile.dat examples/3p02A/domain_list_withPath_reindex  test_out/3p02A 5

### output:

```

<h4> Run SAXSDom-abinitio </h4>

```
Usage:
sh bin/run_SAXSDom_abinitio.sh <target id> <fasta file> <domain list with path> <output directory> <number of cores to run, default:1>

Example:
sh bin/run_SAXSDom_abinitio.sh 3p02A  examples/3p02A/3p02A.fasta  examples/3p02A/domain_list_withPath_reindex  test_out/3p02A_abinitio 5

### output:

```

<h4> Run Modeller </h4>

<h4> Run AIDA </h4>
