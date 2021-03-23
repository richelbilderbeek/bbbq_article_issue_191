# bbbq_article_issues_191

BBBQ article issue 191

## Overview

This is the first sequence and its match to the epitope `MEKSSLTQHSW`.

Scroll to the right to see the match in the sequence,
or see the zoomed-in part below.

```
> TE2IP_HUMAN Telomeric repeat-binding factor 2-interacting protein 1 OS=Homo sapiens (Human) OX=9606 GN=TERF2IP PE=1 SV=1
MAEAMDLGKDPNGPTHSSTLFVRDDGSSMSFYVRPSPAKRRLSTLILHGGGTVCRVQEPGAVLLAQPGEALAEASGDFISTQYILDCVERNERLELEAYRLGPASAADTGSEAKPGALAEGAAEPEPQRHAGRIAFTDADDVAILTYVKENARSPSSVTGNALWKAMEKSSLTQHSWQSLKDRYLKHLRGQEHKYLLGDAPVSPSSQKLKRKAEEDPEAADSGEPQNKRTPDLPEEEYVKEEIQENEEAVKKMLVEATREFEEVVVDESPPDFEIHITMCDDDPPTPEEDSETQPDEEEEEEEEKVSQPEVGAAIKIIRQLMEKFNLDLSTVTQAFLKNSGELEATSAFLASGQRADGYPIWSRQDDIDLQKDDEDTREALVKKFGAQNVARRIEFRKK
                                                                                                                                                            ^        ^^         ^^        ^
                                                                                                                                                            | before ||         ||        |
                                                                                                                                                            + before +|         ||        |
                                                                                                                                                                      + epitope +|        |
                                                                                                                                                                                 + after  +
```

This is the zoomed-in part:

```
SVTGNALWKAMEKSSLTQHSWQSLKDRYLKH
0987654321-----------1234567890
^^^^^^^^^^           ^^^^^^^^^^
||||||||||           ||||||||||
||||||||||           |||||||||+- ten
||||||||||           ||||||||+-- nine
||||||||||           |||||||+--- eight
||||||||||           ||||||+---- seven
||||||||||           |||||+----- six
||||||||||           ||||+------ five
||||||||||           |||+------- four
||||||||||           ||+-------- three
||||||||||           |+--------- two
||||||||||           +---------- one
|||||||||+- minus_one      
||||||||+-- minus_two     
|||||||+--- minus_three    
||||||+---- minus_four    
|||||+----- minus_fix      
||||+------ minus_six     
|||+------- minus_seven     
||+-------- minus_eight    
|+--------- minus_nine      
+---------- minus_ten      
```

## ten_before_epitopes.txt

The ten AAs before the epitope, e.g.:

```
SVTGNALWKAMEKSSLTQHSW
1234567890^         ^
          |         |
          + epitope +
```

## ten_before_epitopes.csv

The counts of the AAs at the ten spots before the epitopes:

```
SVTGNALWKAMEKSSLTQHSW
0987654321-----------
^^^^^^^^^^           ^^^^^^^^^^
||||||||||           ||||||||||
|||||||||+- minus_one      
||||||||+-- minus_two     
|||||||+--- minus_three    
||||||+---- minus_four    
|||||+----- minus_fix      
||||+------ minus_six     
|||+------- minus_seven     
||+-------- minus_eight    
|+--------- minus_nine      
+---------- minus_ten      
```

## ten_after_epitopes.txt

The ten AAs after the epitope, e.g.:

```
MEKSSLTQHSWQSLKDRYLKH
^         ^0987654321
|         |
+ epitope +
```

## ten_after_epitopes.csv

The counts of the AAs at the ten spots after the epitopes:

```
MEKSSLTQHSWQSLKDRYLKH
0987654321-----------1234567890
           ^^^^^^^^^^
           ||||||||||
           |||||||||+- ten
           ||||||||+-- nine
           |||||||+--- eight
           ||||||+---- seven
           |||||+----- six
           ||||+------ five
           |||+------- four
           ||+-------- three
           |+--------- two
           +---------- one
```

