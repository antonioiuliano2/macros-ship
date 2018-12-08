#General FLUKA readme (Last change 7 December 2018)
Set environment variables at the start of each session:  

`export FLUPRO=/home/antonio/FLUKA`  
`export FLUFOR=gfortran`  

To do it from bash:  

`source enviromentfluka.sh`  

Run FLUKA (project MUST be placed in a different folder from FLUKA):
`$FLUPRO/flutil/rfluka -N0 -M1 nomefile.inp`
where N is the start loop and M the final loop. Defaults: N=0 and M=5. We can omitt .inp estension

How to change an user routine: after modification, compile with fff macro:  

$FLUPRO/flutil/fff mgdraw.f  
then link to FLUKA in order to change and executable with lfluka:  
$FLUPRO/flutil/lfluka -m fluka -o miofluka mgdraw.o
Finally launch the executable:  
$FLUPRO/flutil/rfluka -e miofluka -N0 -M1 nomefile.inp  