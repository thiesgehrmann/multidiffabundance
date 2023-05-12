# mda.maaslin2 works

    Code
      res
    Output
      $res
      # A tibble: 144 x 12
         taxa      varia~1 effec~2     se stat   pvalue qvalu~3  qvalue formula method
         <chr>     <chr>     <dbl>  <dbl> <lgl>   <dbl>   <dbl>   <dbl> <chr>   <chr> 
       1 b00d4499~ Report~   1.02  0.0110 NA    2.71e-9 1.95e-7 3.90e-7 ~Repor~ maasl~
       2 bb370945~ Report~  -1.38  0.190  NA    3.54e-4 1.27e-2 2.55e-2 ~Repor~ maasl~
       3 c7ec5769~ Report~   1.85  0.319  NA    1.17e-3 2.82e-2 5.64e-2 ~Repor~ maasl~
       4 f16a387e~ Report~   2.98  0.720  NA    9.02e-3 1.26e-1 2.51e-1 ~Repor~ maasl~
       5 a1b971c0~ Report~   2.35  0.611  NA    8.43e-3 1.26e-1 2.51e-1 ~Repor~ maasl~
       6 ee0fb26f~ Report~  -2.61  0.655  NA    1.05e-2 1.26e-1 2.51e-1 ~Repor~ maasl~
       7 f59a4d48~ Report~  -2.60  0.683  NA    1.27e-2 1.30e-1 2.61e-1 ~Repor~ maasl~
       8 cba5f5b1~ Report~  -2.06  0.648  NA    1.90e-2 1.71e-1 3.43e-1 ~Repor~ maasl~
       9 d86ef5d6~ Report~   1.33  0.477  NA    3.80e-2 3.04e-1 5.54e-1 ~Repor~ maasl~
      10 fcd4f95c~ Report~   0.500 0.224  NA    7.56e-2 3.66e-1 6.10e-1 ~Repor~ maasl~
      # ... with 134 more rows, 2 more variables: n <int>, freq <chr>, and
      #   abbreviated variable names 1: variable, 2: effectsize,
      #   3: qvalue.withinformula
      
      $res.full
      # A tibble: 543 x 10
         taxa          varia~1 effec~2     se  pvalue formula method     n freq  stat 
         <chr>         <chr>     <dbl>  <dbl>   <dbl> <chr>   <chr>  <int> <chr> <lgl>
       1 b00d44992702~ Report~    1.02 0.0110 2.71e-9 ~Repor~ maasl~     8 No:6~ NA   
       2 X91976489627~ Report~    1.17 0.0761 2.12e-5 ~Repor~ maasl~     8 No:6~ NA   
       3 X3677e15d866~ Report~    1.24 0.105  7.97e-5 ~Repor~ maasl~     8 No:6~ NA   
       4 bb370945a677~ Report~   -1.38 0.190  3.54e-4 ~Repor~ maasl~     8 No:6~ NA   
       5 X24dfe6a7325~ Report~    1.43 0.192  6.90e-4 ~Repor~ maasl~     8 No:6~ NA   
       6 c7ec57695a6d~ Report~    1.85 0.319  1.17e-3 ~Repor~ maasl~     8 No:6~ NA   
       7 f16a387e898d~ Report~    2.98 0.720  9.02e-3 ~Repor~ maasl~     8 No:6~ NA   
       8 a1b971c01ce0~ Report~    2.35 0.611  8.43e-3 ~Repor~ maasl~     8 No:6~ NA   
       9 ee0fb26fca3a~ Report~   -2.61 0.655  1.05e-2 ~Repor~ maasl~     8 No:6~ NA   
      10 f59a4d48a03c~ Report~   -2.60 0.683  1.27e-2 ~Repor~ maasl~     8 No:6~ NA   
      # ... with 533 more rows, and abbreviated variable names 1: variable,
      #   2: effectsize
      
      $summary
      $summary$nsig
      # A tibble: 72 x 3
         taxa                             ReportedAntibioticUsage DaysSinceExperimen~1
         <chr>                                              <int>                <int>
       1 b00d44992702bd3743d7f353638d42f9                       1                    0
       2 bb370945a6777f712cfd963c55f2ff54                       1                    0
       3 c7ec57695a6d83b461eba523752755ec                       0                    0
       4 f16a387e898dbae6432b2ed0241120bd                       0                    0
       5 a1b971c01ce0b220275c218adbe717c1                       0                    0
       6 ee0fb26fca3afb05e0b288d6fcab899e                       0                    0
       7 f59a4d48a03cc224462bf83e4d4ab126                       0                    0
       8 cba5f5b1153235425628bd77a28257c7                       0                    0
       9 d86ef5d6394f5dbeb945f39aa25e7426                       0                    0
      10 fcd4f95c05b868060121ff709085bf21                       0                    0
      # ... with 62 more rows, and abbreviated variable name
      #   1: DaysSinceExperimentStart
      
      $summary$pvalue
      # A tibble: 72 x 3
         taxa                             ReportedAntibioticUsage DaysSinceExperimen~1
         <chr>                            <list>                  <list>              
       1 b00d44992702bd3743d7f353638d42f9 <dbl [1]>               <dbl [1]>           
       2 bb370945a6777f712cfd963c55f2ff54 <dbl [1]>               <dbl [1]>           
       3 c7ec57695a6d83b461eba523752755ec <dbl [1]>               <dbl [1]>           
       4 f16a387e898dbae6432b2ed0241120bd <dbl [1]>               <dbl [1]>           
       5 a1b971c01ce0b220275c218adbe717c1 <dbl [1]>               <dbl [1]>           
       6 ee0fb26fca3afb05e0b288d6fcab899e <dbl [1]>               <dbl [1]>           
       7 f59a4d48a03cc224462bf83e4d4ab126 <dbl [1]>               <dbl [1]>           
       8 cba5f5b1153235425628bd77a28257c7 <dbl [1]>               <dbl [1]>           
       9 d86ef5d6394f5dbeb945f39aa25e7426 <dbl [1]>               <dbl [1]>           
      10 fcd4f95c05b868060121ff709085bf21 <dbl [1]>               <dbl [1]>           
      # ... with 62 more rows, and abbreviated variable name
      #   1: DaysSinceExperimentStart
      
      $summary$qvalue
      # A tibble: 72 x 3
         taxa                             ReportedAntibioticUsage DaysSinceExperimen~1
         <chr>                            <list>                  <list>              
       1 b00d44992702bd3743d7f353638d42f9 <dbl [1]>               <dbl [1]>           
       2 bb370945a6777f712cfd963c55f2ff54 <dbl [1]>               <dbl [1]>           
       3 c7ec57695a6d83b461eba523752755ec <dbl [1]>               <dbl [1]>           
       4 f16a387e898dbae6432b2ed0241120bd <dbl [1]>               <dbl [1]>           
       5 a1b971c01ce0b220275c218adbe717c1 <dbl [1]>               <dbl [1]>           
       6 ee0fb26fca3afb05e0b288d6fcab899e <dbl [1]>               <dbl [1]>           
       7 f59a4d48a03cc224462bf83e4d4ab126 <dbl [1]>               <dbl [1]>           
       8 cba5f5b1153235425628bd77a28257c7 <dbl [1]>               <dbl [1]>           
       9 d86ef5d6394f5dbeb945f39aa25e7426 <dbl [1]>               <dbl [1]>           
      10 fcd4f95c05b868060121ff709085bf21 <dbl [1]>               <dbl [1]>           
      # ... with 62 more rows, and abbreviated variable name
      #   1: DaysSinceExperimentStart
      
      $summary$effectsize
      # A tibble: 72 x 3
         taxa                             ReportedAntibioticUsage DaysSinceExperimen~1
         <chr>                            <list>                  <list>              
       1 b00d44992702bd3743d7f353638d42f9 <dbl [1]>               <dbl [1]>           
       2 bb370945a6777f712cfd963c55f2ff54 <dbl [1]>               <dbl [1]>           
       3 c7ec57695a6d83b461eba523752755ec <dbl [1]>               <dbl [1]>           
       4 f16a387e898dbae6432b2ed0241120bd <dbl [1]>               <dbl [1]>           
       5 a1b971c01ce0b220275c218adbe717c1 <dbl [1]>               <dbl [1]>           
       6 ee0fb26fca3afb05e0b288d6fcab899e <dbl [1]>               <dbl [1]>           
       7 f59a4d48a03cc224462bf83e4d4ab126 <dbl [1]>               <dbl [1]>           
       8 cba5f5b1153235425628bd77a28257c7 <dbl [1]>               <dbl [1]>           
       9 d86ef5d6394f5dbeb945f39aa25e7426 <dbl [1]>               <dbl [1]>           
      10 fcd4f95c05b868060121ff709085bf21 <dbl [1]>               <dbl [1]>           
      # ... with 62 more rows, and abbreviated variable name
      #   1: DaysSinceExperimentStart
      
      $summary$method_order
      # A tibble: 1 x 1
        method  
        <chr>   
      1 maaslin2
      
      

