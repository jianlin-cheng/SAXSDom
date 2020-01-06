# SAXSDom
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

**Example**
```
cd examples

sh T1-run-3p02A.sh
```



<h4> Run SAXSDom </h4>

<h4> Run SAXSDom-abinitio </h4>

<h4> Run Modeller </h4>

<h4> Run AIDA </h4>
