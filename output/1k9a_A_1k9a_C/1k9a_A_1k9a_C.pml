reinitialize
load 1k9a_A_1k9a_C.pdb
bg_color white
color grey
select region0, resi 0-72
select region0, region0 + resi 168-195
select region0, region0 + resi 202-215
select region0, region0 + resi 221-433
set_color colour0 = [0  ,0  ,255]
color colour0, region0
select region1, resi 73-165
select region1, region1 + resi 216-217
set_color colour1 = [255,0  ,0  ]
color colour1, region1
select region2, resi 72-73
set_color colour2 = [0  ,255,0  ]
color colour2, region2
select region3, resi 165-169
set_color colour3 = [0  ,255,0  ]
color colour3, region3
select region4, resi 192-216
set_color colour4 = [0  ,255,0  ]
color colour4, region4
select region5, resi 217-222
set_color colour5 = [0  ,255,0  ]
color colour5, region5
set dash_gap, 0
set dash_radius, 0.2
