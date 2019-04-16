#!/bin/bash
awk '
 function bias(arg,laglabel){
 printf("%s_bias: MATHEVAL ARG=%s,%s VAR=x,y FUNC=-x*y*2.5 PERIODIC=NO\n",arg,arg,laglabel);
 printf("BIASVALUE ARG=%s_bias\n",arg);
 }

 function bias_reweight(torsion,kbt,lagmult,type,mult){
  if(mult==1) mult_="" ;else mult_=mult;
  printf("%sbias%s%s: MATHEVAL ARG=kbt,%s,%ssum%s%s FUNC=-x*y*z PERIODIC=NO\n",mult_,type,torsion,lagmult,mult_,type,torsion);
  printf("BIASVALUE ARG=%sbias%s%s\n",mult_,type,torsion);
 }
 
 function bias_reweight_chi(torsion,kbt,lagmult,type,mult){
  if(mult==1) mult_="" ;else mult_=mult;
  printf("%sbias%s%schi: MATHEVAL ARG=kbt,%s,%s%ssum%schi FUNC=-x*y*z PERIODIC=NO\n",mult_,type,torsion,lagmult,mult_,torsion,type);
  printf("BIASVALUE ARG=%sbias%s%schi\n",mult_,type,torsion);
 }
 function strip(s){
  return substr(s, 1, length(s)-1);
 }
BEGIN{
  
  print "# vim: ft=plumed";

  #residues in pdb start from?:
  first_resn=1;
  
  #residues names (in order to choose correct Karplus parameters)
  seq[1]="A";
  seq[2]="A";
  seq[3]="A";
  seq[4]="A";
  #number of residues
  nres=4;

  tmp="";
  for(i=first_resn+1;i<=nres;i++){
    ##alpha penalty 
    printf("alpha%d: TORSION ATOMS=@alpha-%d\n",i,i)
  
#cosalpha1
    arg_cosalpha1=sprintf("%salpha%d,",     arg_cosalpha1,i);
    var_cosalpha1=sprintf("%sv%d,",         var_cosalpha1,i);
   func_cosalpha1=sprintf("%scos(%d*v%d)+",func_cosalpha1,1,i); 
#cosalpha2
    arg_cosalpha2=sprintf("%salpha%d,",     arg_cosalpha2,i);
    var_cosalpha2=sprintf("%sv%d,",         var_cosalpha2,i);
   func_cosalpha2=sprintf("%scos(%d*v%d)+",func_cosalpha2,2,i); 
#cosalpha3
    arg_cosalpha3=sprintf("%salpha%d,",     arg_cosalpha3,i);
    var_cosalpha3=sprintf("%sv%d,",         var_cosalpha3,i);
   func_cosalpha3=sprintf("%scos(%d*v%d)+",func_cosalpha3,3,i); 

#sinalpha1
    arg_sinalpha1=sprintf("%salpha%d,",     arg_sinalpha1,i);
    var_sinalpha1=sprintf("%sv%d,",         var_sinalpha1,i);
   func_sinalpha1=sprintf("%ssin(%d*v%d)+",func_sinalpha1,1,i); 

#sinalpha2
    arg_sinalpha2=sprintf("%salpha%d,",     arg_sinalpha2,i);
    var_sinalpha2=sprintf("%sv%d,",         var_sinalpha2,i);
   func_sinalpha2=sprintf("%ssin(%d*v%d)+",func_sinalpha2,2,i); 
#sinalpha3
    arg_sinalpha3=sprintf("%salpha%d,",     arg_sinalpha3,i);
    var_sinalpha3=sprintf("%sv%d,",         var_sinalpha3,i);
   func_sinalpha3=sprintf("%ssin(%d*v%d)+",func_sinalpha3,3,i); 
#####

##beta
printf("beta%d: TORSION ATOMS=@beta-%d\n",i,i)
#cosbeta1
    arg_cosbeta1=sprintf("%sbeta%d,",     arg_cosbeta1,i);
    var_cosbeta1=sprintf("%sv%d,",         var_cosbeta1,i);
   func_cosbeta1=sprintf("%scos(%d*v%d)+",func_cosbeta1,1,i); 
#cosbeta2
    arg_cosbeta2=sprintf("%sbeta%d,",     arg_cosbeta2,i);
    var_cosbeta2=sprintf("%sv%d,",         var_cosbeta2,i);
   func_cosbeta2=sprintf("%scos(%d*v%d)+",func_cosbeta2,2,i); 
#cosbeta3
    arg_cosbeta3=sprintf("%sbeta%d,",     arg_cosbeta3,i);
    var_cosbeta3=sprintf("%sv%d,",         var_cosbeta3,i);
   func_cosbeta3=sprintf("%scos(%d*v%d)+",func_cosbeta3,3,i); 

#sinbeta1
    arg_sinbeta1=sprintf("%sbeta%d,",     arg_sinbeta1,i);
    var_sinbeta1=sprintf("%sv%d,",         var_sinbeta1,i);
   func_sinbeta1=sprintf("%ssin(%d*v%d)+",func_sinbeta1,1,i); 

#sinbeta2
    arg_sinbeta2=sprintf("%sbeta%d,",     arg_sinbeta2,i);
    var_sinbeta2=sprintf("%sv%d,",         var_sinbeta2,i);
   func_sinbeta2=sprintf("%ssin(%d*v%d)+",func_sinbeta2,2,i); 
#sinbeta3
    arg_sinbeta3=sprintf("%sbeta%d,",     arg_sinbeta3,i);
    var_sinbeta3=sprintf("%sv%d,",         var_sinbeta3,i);
   func_sinbeta3=sprintf("%ssin(%d*v%d)+",func_sinbeta3,3,i); 
#####
   }
   
   ##chi purines
   for(i=first_resn;i<=nres;i++){
    
    if(seq[i]=="A" || seq[i]=="G"){
      printf("chi%d: TORSION ATOMS=@chi-%d\n",i,i)
      #cospurcos1
          arg_cospur1=sprintf("%schi%d,",       arg_cospur1,i);
          var_cospur1=sprintf("%sv%d,",         var_cospur1,i);
         func_cospur1=sprintf("%scos(%d*v%d)+",func_cospur1,1,i); 
      #cospurcos2
          arg_cospur2=sprintf("%schi%d,",       arg_cospur2,i);
          var_cospur2=sprintf("%sv%d,",         var_cospur2,i);
         func_cospur2=sprintf("%scos(%d*v%d)+",func_cospur2,2,i); 
      #cospurcos3
          arg_cospur3=sprintf("%schi%d,",       arg_cospur3,i);
          var_cospur3=sprintf("%sv%d,",         var_cospur3,i);
         func_cospur3=sprintf("%scos(%d*v%d)+",func_cospur3,3,i); 
      
      #sinpurcos1
          arg_sinpur1=sprintf("%schi%d,",       arg_sinpur1,i);
          var_sinpur1=sprintf("%sv%d,",         var_sinpur1,i);
         func_sinpur1=sprintf("%ssin(%d*v%d)+",func_sinpur1,1,i); 
      
      #sinpurcos2
          arg_sinpur2=sprintf("%schi%d,",       arg_sinpur2,i);
          var_sinpur2=sprintf("%sv%d,",         var_sinpur2,i);
         func_sinpur2=sprintf("%ssin(%d*v%d)+",func_sinpur2,2,i); 
      #sinpurcos3
          arg_sinpur3=sprintf("%schi%d,",       arg_sinpur3,i);
          var_sinpur3=sprintf("%sv%d,",         var_sinpur3,i);
         func_sinpur3=sprintf("%ssin(%d*v%d)+",func_sinpur3,3,i); 
      #####
    }else{

#cospurcos1
    arg_cospur1=sprintf("%schi%d,",       arg_cospur1,i);
    var_cospur1=sprintf("%sv%d,",         var_cospur1,i);
   func_cospur1=sprintf("%s0*(v%d)+",    func_cospur1,i); 
#cospurcos2
    arg_cospur2=sprintf("%schi%d,",       arg_cospur2,i);
    var_cospur2=sprintf("%sv%d,",         var_cospur2,i);
   func_cospur2=sprintf("%s0*(v%d)+",    func_cospur2,i); 
#cospurcos3
    arg_cospur3=sprintf("%schi%d,",       arg_cospur3,i);
    var_cospur3=sprintf("%sv%d,",         var_cospur3,i);
   func_cospur3=sprintf("%s0*(v%d)+",    func_cospur3,i); 

#sinpurcos1
    arg_sinpur1=sprintf("%schi%d,",       arg_sinpur1,i);
    var_sinpur1=sprintf("%sv%d,",         var_sinpur1,i);
   func_sinpur1=sprintf("%s0*(v%d)+",    func_sinpur1,i); 

#sinpurcos2
    arg_sinpur2=sprintf("%schi%d,",       arg_sinpur2,i);
    var_sinpur2=sprintf("%sv%d,",         var_sinpur2,i);
   func_sinpur2=sprintf("%s0*(v%d)+",    func_sinpur2,i); 
#sinpurcos3
    arg_sinpur3=sprintf("%schi%d,",       arg_sinpur3,i);
    var_sinpur3=sprintf("%sv%d,",         var_sinpur3,i);
   func_sinpur3=sprintf("%s0*(v%d)+",    func_sinpur3,i); 
#####

}
    if(seq[i]=="C" || seq[i]=="U"){

      printf("chi%d: TORSION ATOMS=@chi-%d\n",i,i)
      #cospyrcos1
          arg_cospyr1=sprintf("%schi%d,",       arg_cospyr1,i);
          var_cospyr1=sprintf("%sv%d,",         var_cospyr1,i);
         func_cospyr1=sprintf("%scos(%d*v%d)+",func_cospyr1,1,i); 
      #cospyrcos2pyr
          arg_cospyr2=sprintf("%schi%d,",       arg_cospyr2,i);
          var_cospyr2=sprintf("%sv%d,",         var_cospyr2,i);
         func_cospyr2=sprintf("%scos(%d*v%d)+",func_cospyr2,2,i); 
      #cospyrcos3pyr
          arg_cospyr3=sprintf("%schi%d,",       arg_cospyr3,i);
          var_cospyr3=sprintf("%sv%d,",         var_cospyr3,i);
         func_cospyr3=sprintf("%scos(%d*v%d)+",func_cospyr3,3,i); 
      
      #sinpyrcos1pyr
          arg_sinpyr1=sprintf("%schi%d,",       arg_sinpyr1,i);
          var_sinpyr1=sprintf("%sv%d,",         var_sinpyr1,i);
         func_sinpyr1=sprintf("%ssin(%d*v%d)+",func_sinpyr1,1,i); 
      
      #sinpyrcos2pyr
          arg_sinpyr2=sprintf("%schi%d,",       arg_sinpyr2,i);
          var_sinpyr2=sprintf("%sv%d,",         var_sinpyr2,i);
         func_sinpyr2=sprintf("%ssin(%d*v%d)+",func_sinpyr2,2,i); 
      #sinpyrcos3pyr
          arg_sinpyr3=sprintf("%schi%d,",       arg_sinpyr3,i);
          var_sinpyr3=sprintf("%sv%d,",         var_sinpyr3,i);
         func_sinpyr3=sprintf("%ssin(%d*v%d)+",func_sinpyr3,3,i); 
      #####
      
    }else{
 #
    arg_cospyr1=sprintf("%schi%d,",       arg_cospyr1,i);
    var_cospyr1=sprintf("%sv%d,",         var_cospyr1,i);
   func_cospyr1=sprintf("%s0*(v%d)+",    func_cospyr1,i); 
#
    arg_cospyr2=sprintf("%schi%d,",       arg_cospyr2,i);
    var_cospyr2=sprintf("%sv%d,",         var_cospyr2,i);
   func_cospyr2=sprintf("%s0*(v%d)+",    func_cospyr2,i); 
#
    arg_cospyr3=sprintf("%schi%d,",       arg_cospyr3,i);
    var_cospyr3=sprintf("%sv%d,",         var_cospyr3,i);
   func_cospyr3=sprintf("%s0*(v%d)+",    func_cospyr3,i); 

#
    arg_sinpyr1=sprintf("%schi%d,",       arg_sinpyr1,i);
    var_sinpyr1=sprintf("%sv%d,",         var_sinpyr1,i);
   func_sinpyr1=sprintf("%s0*(v%d)+",    func_sinpyr1,i); 

#
    arg_sinpyr2=sprintf("%schi%d,",       arg_sinpyr2,i);
    var_sinpyr2=sprintf("%sv%d,",         var_sinpyr2,i);
   func_sinpyr2=sprintf("%s0*(v%d)+",    func_sinpyr2,i); 
#
    arg_sinpyr3=sprintf("%schi%d,",       arg_sinpyr3,i);
    var_sinpyr3=sprintf("%sv%d,",         var_sinpyr3,i);
   func_sinpyr3=sprintf("%s0*(v%d)+",    func_sinpyr3,i); 
#####

    
    }
  } 

for(i=first_resn;i<=nres-1;i++){

    ##zeta penalty 
     printf("zeta%d: TORSION ATOMS=@zeta-%d\n",i,i)
     #coszeta1
    arg_coszeta1=sprintf("%szeta%d,",     arg_coszeta1,i);
    var_coszeta1=sprintf("%sv%d,",         var_coszeta1,i);
   func_coszeta1=sprintf("%scos(%d*v%d)+",func_coszeta1,1,i); 
#coszeta2
    arg_coszeta2=sprintf("%szeta%d,",     arg_coszeta2,i);
    var_coszeta2=sprintf("%sv%d,",         var_coszeta2,i);
   func_coszeta2=sprintf("%scos(%d*v%d)+",func_coszeta2,2,i); 
#coszeta3
    arg_coszeta3=sprintf("%szeta%d,",     arg_coszeta3,i);
    var_coszeta3=sprintf("%sv%d,",         var_coszeta3,i);
   func_coszeta3=sprintf("%scos(%d*v%d)+",func_coszeta3,3,i); 

#sinzeta1
    arg_sinzeta1=sprintf("%szeta%d,",     arg_sinzeta1,i);
    var_sinzeta1=sprintf("%sv%d,",         var_sinzeta1,i);
   func_sinzeta1=sprintf("%ssin(%d*v%d)+",func_sinzeta1,1,i); 

#sinzeta2
    arg_sinzeta2=sprintf("%szeta%d,",     arg_sinzeta2,i);
    var_sinzeta2=sprintf("%sv%d,",         var_sinzeta2,i);
   func_sinzeta2=sprintf("%ssin(%d*v%d)+",func_sinzeta2,2,i); 
#sinzeta3
    arg_sinzeta3=sprintf("%szeta%d,",     arg_sinzeta3,i);
    var_sinzeta3=sprintf("%sv%d,",         var_sinzeta3,i);
   func_sinzeta3=sprintf("%ssin(%d*v%d)+",func_sinzeta3,3,i); 
#####
printf("epsilon%d: TORSION ATOMS=@epsilon-%d\n",i,i)
  #cosepsilon1
    arg_cosepsilon1=sprintf("%sepsilon%d,",     arg_cosepsilon1,i);
    var_cosepsilon1=sprintf("%sv%d,",         var_cosepsilon1,i);
   func_cosepsilon1=sprintf("%scos(%d*v%d)+",func_cosepsilon1,1,i); 
#cosepsilon2
    arg_cosepsilon2=sprintf("%sepsilon%d,",     arg_cosepsilon2,i);
    var_cosepsilon2=sprintf("%sv%d,",         var_cosepsilon2,i);
   func_cosepsilon2=sprintf("%scos(%d*v%d)+",func_cosepsilon2,2,i); 
#cosepsilon3
    arg_cosepsilon3=sprintf("%sepsilon%d,",     arg_cosepsilon3,i);
    var_cosepsilon3=sprintf("%sv%d,",         var_cosepsilon3,i);
   func_cosepsilon3=sprintf("%scos(%d*v%d)+",func_cosepsilon3,3,i); 

#sinepsilon1
    arg_sinepsilon1=sprintf("%sepsilon%d,",     arg_sinepsilon1,i);
    var_sinepsilon1=sprintf("%sv%d,",         var_sinepsilon1,i);
   func_sinepsilon1=sprintf("%ssin(%d*v%d)+",func_sinepsilon1,1,i); 

#sinepsilon2
    arg_sinepsilon2=sprintf("%sepsilon%d,",     arg_sinepsilon2,i);
    var_sinepsilon2=sprintf("%sv%d,",         var_sinepsilon2,i);
   func_sinepsilon2=sprintf("%ssin(%d*v%d)+",func_sinepsilon2,2,i); 
#sinepsilon3
    arg_sinepsilon3=sprintf("%sepsilon%d,",     arg_sinepsilon3,i);
    var_sinepsilon3=sprintf("%sv%d,",         var_sinepsilon3,i);
   func_sinepsilon3=sprintf("%ssin(%d*v%d)+",func_sinepsilon3,3,i); 
#####
 
  }

  for(i=first_resn;i<=nres;i++){
  
    #gamma
    printf("gamma%d: TORSION ATOMS=@gamma-%d\n",i,i)
   #cosgamma1
    arg_cosgamma1=sprintf("%sgamma%d,",     arg_cosgamma1,i);
    var_cosgamma1=sprintf("%sv%d,",         var_cosgamma1,i);
   func_cosgamma1=sprintf("%scos(%d*v%d)+",func_cosgamma1,1,i); 
#cosgamma2
    arg_cosgamma2=sprintf("%sgamma%d,",     arg_cosgamma2,i);
    var_cosgamma2=sprintf("%sv%d,",         var_cosgamma2,i);
   func_cosgamma2=sprintf("%scos(%d*v%d)+",func_cosgamma2,2,i); 
#cosgamma3
    arg_cosgamma3=sprintf("%sgamma%d,",     arg_cosgamma3,i);
    var_cosgamma3=sprintf("%sv%d,",         var_cosgamma3,i);
   func_cosgamma3=sprintf("%scos(%d*v%d)+",func_cosgamma3,3,i); 

#singamma1
    arg_singamma1=sprintf("%sgamma%d,",     arg_singamma1,i);
    var_singamma1=sprintf("%sv%d,",         var_singamma1,i);
   func_singamma1=sprintf("%ssin(%d*v%d)+",func_singamma1,1,i); 

#singamma2
    arg_singamma2=sprintf("%sgamma%d,",     arg_singamma2,i);
    var_singamma2=sprintf("%sv%d,",         var_singamma2,i);
   func_singamma2=sprintf("%ssin(%d*v%d)+",func_singamma2,2,i); 
#singamma3
    arg_singamma3=sprintf("%sgamma%d,",     arg_singamma3,i);
    var_singamma3=sprintf("%sv%d,",         var_singamma3,i);
   func_singamma3=sprintf("%ssin(%d*v%d)+",func_singamma3,3,i); 
#####
 

     #delta
    printf("delta%d: TORSION ATOMS=@delta-%d\n",i,i)
#cosdelta1
    arg_cosdelta1=sprintf("%sdelta%d,",     arg_cosdelta1,i);
    var_cosdelta1=sprintf("%sv%d,",         var_cosdelta1,i);
   func_cosdelta1=sprintf("%scos(%d*v%d)+",func_cosdelta1,1,i); 
#cosdelta2
    arg_cosdelta2=sprintf("%sdelta%d,",     arg_cosdelta2,i);
    var_cosdelta2=sprintf("%sv%d,",         var_cosdelta2,i);
   func_cosdelta2=sprintf("%scos(%d*v%d)+",func_cosdelta2,2,i); 
#cosdelta3
    arg_cosdelta3=sprintf("%sdelta%d,",     arg_cosdelta3,i);
    var_cosdelta3=sprintf("%sv%d,",         var_cosdelta3,i);
   func_cosdelta3=sprintf("%scos(%d*v%d)+",func_cosdelta3,3,i); 

#sindelta1
    arg_sindelta1=sprintf("%sdelta%d,",     arg_sindelta1,i);
    var_sindelta1=sprintf("%sv%d,",         var_sindelta1,i);
   func_sindelta1=sprintf("%ssin(%d*v%d)+",func_sindelta1,1,i); 

#sindelta2
    arg_sindelta2=sprintf("%sdelta%d,",     arg_sindelta2,i);
    var_sindelta2=sprintf("%sv%d,",         var_sindelta2,i);
   func_sindelta2=sprintf("%ssin(%d*v%d)+",func_sindelta2,2,i); 
#sindelta3
    arg_sindelta3=sprintf("%sdelta%d,",     arg_sindelta3,i);
    var_sindelta3=sprintf("%sv%d,",         var_sindelta3,i);
   func_sindelta3=sprintf("%ssin(%d*v%d)+",func_sindelta3,3,i); 
#####


  }
  
  #printf("pursumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospur1),strip(var_cospur1),strip(func_cospur1));
  #printf("2pursumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospur2),strip(var_cospur2),strip(func_cospur2));
  #printf("pyrsumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospyr1),strip(var_cospyr1),strip(func_cospyr1));

printf("sumcosalpha: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosalpha1),strip(var_cosalpha1),strip(func_cosalpha1));
printf("sumsinalpha: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinalpha1),strip(var_sinalpha1),strip(func_sinalpha1));

printf("2sumcosalpha: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosalpha2),strip(var_cosalpha2),strip(func_cosalpha2));
printf("2sumsinalpha: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinalpha2),strip(var_sinalpha2),strip(func_sinalpha2));

printf("3sumcosalpha: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosalpha3),strip(var_cosalpha3),strip(func_cosalpha3));
printf("3sumsinalpha: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinalpha3),strip(var_sinalpha3),strip(func_sinalpha3));

printf("sumcosbeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosbeta1),strip(var_cosbeta1),strip(func_cosbeta1));
printf("sumsinbeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinbeta1),strip(var_sinbeta1),strip(func_sinbeta1));

printf("2sumcosbeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosbeta2),strip(var_cosbeta2),strip(func_cosbeta2));
printf("2sumsinbeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinbeta2),strip(var_sinbeta2),strip(func_sinbeta2));

printf("3sumcosbeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosbeta3),strip(var_cosbeta3),strip(func_cosbeta3));
printf("3sumsinbeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinbeta3),strip(var_sinbeta3),strip(func_sinbeta3));

printf("sumcosgamma: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosgamma1),strip(var_cosgamma1),strip(func_cosgamma1));
printf("sumsingamma: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_singamma1),strip(var_singamma1),strip(func_singamma1));

printf("2sumcosgamma: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosgamma2),strip(var_cosgamma2),strip(func_cosgamma2));
printf("2sumsingamma: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_singamma2),strip(var_singamma2),strip(func_singamma2));

printf("3sumcosgamma: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosgamma3),strip(var_cosgamma3),strip(func_cosgamma3));
printf("3sumsingamma: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_singamma3),strip(var_singamma3),strip(func_singamma3));

printf("sumcosdelta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosdelta1),strip(var_cosdelta1),strip(func_cosdelta1));
printf("sumsindelta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sindelta1),strip(var_sindelta1),strip(func_sindelta1));

printf("2sumcosdelta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosdelta2),strip(var_cosdelta2),strip(func_cosdelta2));
printf("2sumsindelta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sindelta2),strip(var_sindelta2),strip(func_sindelta2));

printf("3sumcosdelta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosdelta3),strip(var_cosdelta3),strip(func_cosdelta3));
printf("3sumsindelta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sindelta3),strip(var_sindelta3),strip(func_sindelta3));

printf("sumcosepsilon: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosepsilon1),strip(var_cosepsilon1),strip(func_cosepsilon1));
printf("sumsinepsilon: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinepsilon1),strip(var_sinepsilon1),strip(func_sinepsilon1));

printf("2sumcosepsilon: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosepsilon2),strip(var_cosepsilon2),strip(func_cosepsilon2));
printf("2sumsinepsilon: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinepsilon2),strip(var_sinepsilon2),strip(func_sinepsilon2));

printf("3sumcosepsilon: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cosepsilon3),strip(var_cosepsilon3),strip(func_cosepsilon3));
printf("3sumsinepsilon: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinepsilon3),strip(var_sinepsilon3),strip(func_sinepsilon3));

printf("sumcoszeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_coszeta1),strip(var_coszeta1),strip(func_coszeta1));
printf("sumsinzeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinzeta1),strip(var_sinzeta1),strip(func_sinzeta1));

printf("2sumcoszeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_coszeta2),strip(var_coszeta2),strip(func_coszeta2));
printf("2sumsinzeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinzeta2),strip(var_sinzeta2),strip(func_sinzeta2));

printf("3sumcoszeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_coszeta3),strip(var_coszeta3),strip(func_coszeta3));
printf("3sumsinzeta: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinzeta3),strip(var_sinzeta3),strip(func_sinzeta3));

 printf("pursumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospur1),strip(var_cospur1),strip(func_cospur1));
printf("2pursumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospur2),strip(var_cospur2),strip(func_cospur2));
printf("3pursumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospur3),strip(var_cospur3),strip(func_cospur3));

 printf("pursumsinchi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinpur1),strip(var_sinpur1),strip(func_sinpur1));
printf("2pursumsinchi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinpur2),strip(var_sinpur2),strip(func_sinpur2));
printf("3pursumsinchi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinpur3),strip(var_sinpur3),strip(func_sinpur3));

 printf("pyrsumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospyr1),strip(var_cospyr1),strip(func_cospyr1));
printf("2pyrsumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospyr2),strip(var_cospyr2),strip(func_cospyr2));
printf("3pyrsumcoschi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_cospyr3),strip(var_cospyr3),strip(func_cospyr3));

 printf("pyrsumsinchi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinpyr1),strip(var_sinpyr1),strip(func_sinpyr1));
printf("2pyrsumsinchi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinpyr2),strip(var_sinpyr2),strip(func_sinpyr2));
printf("3pyrsumsinchi: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",strip(arg_sinpyr3),strip(var_sinpyr3),strip(func_sinpyr3));
 
  add_bias=1
  if(add_bias){

lag[1]=0.0287429
lag[2]=-0.118088
lag[3]=0.36924
lag[4]=-0.0688751
lag[5]=-0.0617793
lag[6]=-0.0835629
lag[7]=0.04596
lag[8]=-0.0857072
lag[9]=0.0806448
lag[10]=-0.0114061
lag[11]=0.0389671
lag[12]=0.220101
lag[13]=-0.0078801
lag[14]=0.0677268
lag[15]=-0.0864179
lag[16]=-0.0144299
lag[17]=-0.135683
lag[18]=0.199287
lag[19]=-0.269494
lag[20]=0.0358081
lag[21]=0.0866801
lag[22]=0.0916538
lag[23]=0.0646414
lag[24]=0.00892007
lag[25]=0.190771
lag[26]=-0.109477
lag[27]=0.0540643
lag[28]=0.150732
lag[29]=-0.0545985
lag[30]=-0.116106
lag[31]=-0.0731506
lag[32]=-0.0460672
lag[33]=0.122662
lag[34]=-0.0649582
lag[35]=-0.00194904
lag[36]=0.10993
lag[37]=0.0424341
lag[38]=0.0616814
lag[39]=0.0321497
lag[40]=0.0160872
lag[41]=0.0209927
lag[42]=0.0510608
lag[43]=-0.0777578
lag[44]=0.0159556
lag[45]=0.0169775
lag[46]=0.0374566
lag[47]=0.143994
lag[48]=0.140334

dim=48;
    printf("kbt: CONSTANT VALUE=2.5\n")
    biasarg="kbt,"
    biasarg=sprintf("%s%s",biasarg,"sumcosalpha,sumcosbeta,sumcosgamma,sumcosdelta,sumcosepsilon,sumcoszeta,sumsinalpha,sumsinbeta,sumsingamma,sumsindelta,sumsinepsilon,sumsinzeta,pursumcoschi,pyrsumcoschi,pursumsinchi,pyrsumsinchi,2sumcosalpha,2sumcosbeta,2sumcosgamma,2sumcosdelta,2sumcosepsilon,2sumcoszeta,2sumsinalpha,2sumsinbeta,2sumsingamma,2sumsindelta,2sumsinepsilon,2sumsinzeta,2pursumcoschi,2pyrsumcoschi,2pursumsinchi,2pyrsumsinchi,3sumcosalpha,3sumcosbeta,3sumcosgamma,3sumcosdelta,3sumcosepsilon,3sumcoszeta,3sumsinalpha,3sumsinbeta,3sumsingamma,3sumsindelta,3sumsinepsilon,3sumsinzeta,3pursumcoschi,3pyrsumcoschi,3pursumsinchi,3pyrsumsinchi")
    biasvar="kbt,"
    biasfunc="kbt*("
    for(i=1;i<=dim;i++){
      biasvar=sprintf("%sv%d,",biasvar,i)
      biasfunc=sprintf("%s(%f*v%d)+",biasfunc,lag[i],i)
    }
    biasvar=strip(biasvar)
    biasfunc=strip(biasfunc)
    biasfunc=sprintf("%s)",biasfunc) 
    printf("mybias: MATHEVAL ARG=%s VAR=%s FUNC=%s PERIODIC=NO\n",biasarg,biasvar,biasfunc);
    printf("BIASVALUE ARG=mybias\n") 
      

  } 

}'

