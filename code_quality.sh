#!/bin/bash
cq_result=$(radon cc --min D -e *pyx genmechanics)
echo $cq_result
if [[ "$cq_result" ]];
  then 
	  echo -e "Error in code quality check, run radon to simplify functions">&2;
	  exit 64;
	
fi;
