# bbbq_article_issues_191

BBBQ article issue 191

## ten_before_epitopes.txt

The ten AAs before the epitope, e.g.:

SVTGNALWKAMEKSSLTQHSW
1234567890^         ^
          |         |
          + epitope +

## ten_before_epitopes.csv

The counts of the AAs at the ten spots before the epitopes:

SVTGNALWKAMEKSSLTQHSW
1234567890-----------
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

## ten_after_epitopes.txt

The ten AAs after the epitope, e.g.:

MEKSSLTQHSWQSLKDRYLKH
^         ^0987654321
|         |
+ epitope +

## ten_after_epitopes.csv

The counts of the AAs at the ten spots after the epitopes:

MEKSSLTQHSWQSLKDRYLKH
ALWKAMEKSSLTQHSW
-----------0987654321
           ^^^^^^^^^^
           ||||||||||
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

