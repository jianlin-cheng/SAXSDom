# SAXSDom
This is a bioinformatics tool to use SAXS data to assemble protein domain structures into full-length structural models


**Installation**
```
perl setup_env.pl

perl configure.pl

sh compile_SAXSDom.sh

```



**Example**
export LD_LIBRARY_PATH=/home/jh7x3/IMP2.6/lib:/home/jh7x3/boost_1_55_0/lib:$LD_LIBRARY_PATH

sh ./run_SAXSDom.sh T0998 input output ' -t   -g test_assembly  -d 1 -x  1  --scoreWeight 30_700_700_700 --scoreWeightInitial 30_700_700_700  --
scoreCombine' 50



<h4> Run SAXSDom </h4>
<h4> Run SAXSDom-abinitio </h4>
<h4> Run Modeller </h4>
<h4> Run AIDA </h4>
