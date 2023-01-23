<h1 style="text-align:center;">Trigger Count for Large LIGO Data </h1>

<ol><li> Trigger.py: Main file with computation
   <ul><li> Change env
       <li> Slow GWOSC strain download using wget </ul>
      
   <li> Data:
   <ol><li> url.csv: GWOSC URLs with full cycle LIGO O3b data
       <li> linS.hdf5: Template files (ids and params in attributes) </ol>
     
   <li> Test:
   <ol> <li> test.sh:
        <ul><li> Testing shell script with 3 strain files (4096 s at 256 Hz) instead of 320
            <li> Change env</ul>
        <li> trig_test.sub: Submit file for test.sh</ol>
   
  <li> Run:
   <ol> <li> run.sh:
        <ul><li> shell script with 320 strain files (4096 s at 256 Hz)
            <li> Change env</ul>
        <li> trig_run.sub: Submit file for run.sh</ol>
        
     
