# Sample PLUMED input file for AAAA tetranucleotide
This sample plumed.dat file can be used to apply the fitted corrections, on-the-fly, durin a MD simulation.
## Apply Corrections during simulations
The plumed.dat file can be used to apply the fitted corrections on-the-fly during a MD simulation using PLUMED.

## Generating the bias file
 - Specify sequence in the file gen_bias.sh by modifying this part:
      ```
      seq[1]="A";
      seq[2]="A";
      seq[3]="A";
      seq[4]="A";
     ```  
     For larger systems add `seq[5]` etc.
 - run `bash gen_bias.sh > bias`
 
## Fitted parameters 
They are specified at the end of the gen_bias.sh script and are identified with `lag[1], lag[2]`, etc.
