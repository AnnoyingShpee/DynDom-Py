import numpy as np
import itertools


# def split_consecutive_groups(lst):
#     return [list(group) for key, group in itertools.groupby(lst)]
#
# # initialize list
# test_list = [1, 4, 4, 5, 5, 5, 7, 7, 8, 8, 8, 5, 5]
#
# # printing original list
# print("The original list is : " + str(test_list))
#
# # Identical Consecutive Grouping in list
# # using itertools.groupby()
# res = split_consecutive_groups(test_list)
#
# # printing result
# print("List after grouping is : " + str(res))


# x=np.array([1, 2, 2, 2, 3, 3, 4, 4, 6, 6, 6, 6, 6, 6, 8, 6, 6])
# ranges=[list(g) for _, g in itertools.groupby(range(len(x)), lambda idx:x[idx])]
# # [[0], [1, 2, 3], [4, 5], [6, 7], [8, 9, 10, 11, 12, 13], [14]]
# #   1    2  2  2    3  3    4  4    6  6  6   6   6   6     8
#
# final=[[r[0],r[-1]] for r in ranges if len(r)>0]
# print(final)
# print(type(final))
# print(np.array(final))
# print(type(np.array(final)))
# # [[1, 3], [4, 5], [6, 7], [8, 13]]

# # # Define the matrix B
# B = np.array([
#     [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
#     [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
# ])
#
# R = np.array([
#     [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
# ])
#
# T = np.array([
#     [1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
#     [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
#     [1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
#     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
# ])
#
#
# # Perform row reduction
# def row_reduction(matrix):
#     rows = matrix.shape[0]
#     connections_list: list = []
#     # used_indices = []
#     # Go through each row sequentially starting from the row with the highest index to the lowest index
#     for m in range(rows-1, -1, -1):
#         # print(f"=================================")
#         # print(f"m = {m}")
#         # print(f"---------------------------------")
#         connected_set = set()
#         lowest_index = m
#         or_wise_result = matrix[m]
#         # Go through each column sequentially starting from the current row number to the lowest index
#         for n in range(m, -1, -1):
#             if matrix[m][n] == 1:
#                 list_to_or_wise = matrix[n]
#                 # Perform OR wise logical operation on the 2 arrays
#                 or_wise_result |= list_to_or_wise
#                 print(f"Or wise = {or_wise_result}")
#                 # Returns a tuple (a, b) where a is the list of indices of or_wise_result where the element is 1.
#                 # b is empty.
#                 results = np.where(or_wise_result == 1)
#                 indices_array = results[0]
#                 print(f"Indices = {indices_array}")
#                 # connected_set.update(set(indices_array))
#                 # print(f"Connected set = {connected_set}")
#                 print(f"Connected list = {connections_list}")
#                 ind = [i for i in range(len(connections_list)) if any(np.isin(indices_array, connections_list[i]))]
#                 print(f"ind = {ind}")
#                 if len(ind) > 0:
#                     temp = np.append(connections_list[ind[0]], indices_array)
#                     connections_list[ind[0]] = np.unique(temp, return_counts=False)
#                     # connections_list[ind[0]] = np.append(connections_list[ind[0]], indices_array)
#                     pass
#                 else:
#                     connections_list.append(indices_array)
#
#     return connections_list
#
#
# # Perform row reduction on matrix B
# print(row_reduction(B))

# A = np.empty(shape=[0, 2], dtype=int)
# print(f"Before = {A}")
# V = np.append(A, [3, 4], axis=0)
# print(f"After = {V}")

# A = np.array([[]], dtype=int)
# print(f"Before = {A}")
# print(A.shape)
# A = np.append(A, [[3, 4]], axis=0 if A.shape[1] > 0 else 1)
# # A = np.append(A, [[3, 4]], axis=1)
# print(f"After = {A}")
# print(A.shape)
# # A = np.append(A, [[5, 6]], axis=0)
# A = np.append(A, [[5, 6]], axis=0 if A.shape[1] > 0 else 1)
# print(f"After = {A}")
# print(A.shape)


# A = np.array([[5, 6], [2, 3]], dtype=int)
# print(f"SUM = {S}")





# Generate sample data
# X = np.array([[-58.449466568941546, -50.18029757464576, 84.73776302222632],
# [-56.89368387916321, -62.89381346682054, 63.703058267840696],
# [-116.13827637230091, -61.65467596487683, 107.67054734229158],
# [115.57163346585139, 61.847965927147655, -85.65270339699856],
# [126.66174718454631, 56.73097738970182, -76.47315272118536],
# [130.31496574111455, 51.86764772860048, -76.1222569411798],
# [130.4355744038591, 49.52885849220227, -75.83360151565459],
# [129.6204934482859, 49.739947023525424, -74.74036275912323],
# [130.16910938514104, 49.722179947904074, -74.99278562053023],
# [127.38962674162173, 49.86023016065446, -80.37666612416106],
# [126.38989231058196, 47.729220928829726, -80.4885754054633],
# [121.98942338932878, 46.51449967437468, -78.48927923276825],
# [121.0319575830847, 47.69388809474983, -76.34295251967812],
# [121.04670906839357, 49.14298912133817, -75.10501268487221],
# [121.77937164559278, 52.31425940267442, -78.33719096421721],
# [121.99435074274058, 53.938212769361506, -80.40067054193959],
# [119.26232795810836, 54.139788509080056, -83.78685509564487],
# [118.78697027366186, 50.65336673394709, -85.24669237648725],
# [117.82412504689486, 48.09027007331316, -83.41395938866576],
# [117.20375297967591, 48.55590705230655, -84.1025737787695],
# [118.30192312206523, 49.280085345136975, -82.14097114380316],
# [118.26080080207646, 49.10159201443881, -82.27336756359674],
# [120.99876635106725, 49.08765459774292, -85.50124770906321],
# [128.2000388557561, 48.37950334381356, -87.3449642989589],
# [128.89675754222537, 54.39613034127019, -83.21240849282518],
# [121.76508274050819, 69.31641497697149, -83.37201565474336],
# [118.41708804161723, 70.01612613919599, -83.34928181506473],
# [113.64007908633386, 66.15445895633371, -81.483629819029],
# [113.8403902046481, 51.98183543116537, -86.12713588085147],
# [121.37361020173707, 42.607986749716204, -88.13051881246707],
# [121.17058923108006, 43.92096604340657, -84.72312414763542],
# [122.3756631408534, 49.24887952515407, -83.81058594850202],
# [121.59873958459093, 49.46313161778007, -81.8304361288412],
# [119.90885055232195, 48.37402567184977, -83.88843198789459],
# [117.50182829612424, 50.386168161147644, -89.52747045820207],
# [115.45364381978003, 51.97648031460279, -90.9215050351323],
# [115.76436671561041, 52.29502190424573, -90.69321511734468],
# [114.28311770657574, 53.50594570947649, -88.6904615317605],
# [115.6462148664363, 52.97358190704073, -85.20910695102224],
# [119.23032915506478, 49.52808601246662, -82.83589822932184],
# [122.8140012953622, 51.47909326555544, -83.28048568705836],
# [122.57345676634242, 57.443144815523894, -85.95557158686083],
# [124.10635241273653, 51.158870259118075, -86.02553519581421],
# [124.5973519986448, 49.868107067738514, -83.68079698089777],
# [125.68360217448476, 45.05521612562637, -81.65137939332391],
# [123.11051999929946, 41.9593397625631, -83.67780162751525],
# [121.27190390430414, 45.04656551359326, -84.60533979821533],
# [118.82283185106003, 47.11533361709617, -88.25188493951659],
# [117.66613158279714, 45.88838818443697, -89.52419248315739],
# [124.2568491567966, 53.320901142768676, -87.95752698295804],
# [123.15284983046377, 53.181372431381334, -84.10114655164985],
# [123.41189766092764, 50.571716404145455, -81.82471934419364],
# [127.89557077601081, 52.13234966834012, -80.97093477073244],
# [133.5008553447141, 43.459325370069536, -80.18169331671743],
# [137.76762544370732, 38.46770677874542, -78.27764958702613],
# [143.97217408969564, 43.27224990869659, -73.97950705676114],
# [145.73679694339933, 45.71275849806165, -73.00707549893255],
# [148.35734288817042, 35.50662279265898, -83.8982203372949],
# [149.34435745906967, 37.15753330579319, -88.53501484163033],
# [147.16728139757504, 37.70715266609796, -84.08246513883846],
# [144.47828761182933, 38.044626856109694, -79.40552249699134],
# [136.99025623127542, 45.343959798763606, -72.22070671594956],
# [133.12441960016224, 45.62013909025202, -69.34748014658913],
# [128.6628677233914, 50.50645538705939, -70.67948741013662],
# [128.23858530835238, 54.281327704037544, -66.85997947003077],
# [127.2274005295475, 53.074158484399874, -72.3401079633886],
# [134.28989791071007, 44.651335356791456, -76.62275894527117],
# [122.67268686226332, 50.41998249781646, -80.98799735829579],
# [125.93639437476229, 46.334282937003564, -84.94033619521873],
# [127.17072813923191, 50.934296898267434, -84.95232699012269],
# [127.81033146085056, 47.14164347382415, -82.60551899571274],
# [128.48388287898112, 46.58233730668657, -81.92449684543786],
# [129.54422222018277, 46.47360069998939, -78.49123015393172],
# [130.40728525583165, 46.63758554241213, -77.13357942717934],
# [128.44563981441854, 53.197067984771884, -75.90792442224416],
# [129.50436341982515, 55.16656854775541, -76.50249761966971],
# [127.58788698471356, 52.998135903918424, -77.9182520163073],
# [127.10055865799313, 51.297625901317154, -77.27446285330875],
# [128.31070285215893, 48.45372036347151, -75.94741321090332],
# [127.21608870484026, 47.78574461485217, -76.02530048065032],
# [127.36682929653512, 50.12806374627676, -77.97004700879758],
# [127.813515472154, 52.0646189900577, -77.26276062489242],
# [127.97229920060447, 53.22479291601453, -78.12527084098842],
# [127.0309602200657, 51.58481608887475, -77.63552914331393],
# [126.78602787847942, 50.84003056796237, -77.52386090940911],
# [127.97321696712876, 48.34787218728258, -80.05503077200909],
# [127.69896640294417, 46.03708658820346, -81.92763952680238],
# [124.48004255526673, 49.958042017071406, -83.33216328024065],
# [126.28387020965434, 49.374432787179146, -80.13714661792604],
# [126.62362665548311, 50.23436628955607, -77.49336108992823],
# [126.38297889073752, 49.19418233724262, -76.4843988126941],
# [126.8136010165716, 46.99558384023988, -78.18724890588875],
# [123.42974247074848, 44.19284454732669, -82.64018792274935],
# [122.77986616423888, 44.882772730740335, -82.02707932101728],
# [123.09471688253075, 48.37319062285816, -82.41867210964136],
# [124.33704158909899, 50.80778638127303, -79.89594740508973],
# [126.87269769380212, 48.194325261469174, -79.38707424734643],
# [129.00987518680367, 46.6187809936035, -80.42037795957546],
# [125.38288536959399, 45.534616435497746, -80.34928609805682],
# [125.4513432260803, 47.176100226412466, -79.78355814610686],
# [129.29645434809837, 50.38392875675622, -77.72158165987412],
# [128.08951705188224, 47.964555936672255, -77.23191175484946],
# [129.48914944768458, 46.56230022091844, -77.08871602084442],
# [128.83534246722436, 49.50010937340214, -76.85505834280185],
# [127.00279906558481, 52.348368673354635, -78.66643255080395],
# [121.99391397572036, 54.4043304545715, -79.8557893734274],
# [119.45238496507852, 53.03805896525692, -78.7157706281261],
# [120.09754327472888, 50.21966467755875, -76.96043087815012],
# [122.89175152617409, 46.65709668519527, -74.38263961372007],
# [125.8191858916484, 46.648296595840456, -73.76568934456387],
# [127.50477112897171, 46.232272656915335, -75.19047842397501],
# [128.06788033624082, 47.0193479731321, -76.22868834649287],
# [125.49526974781962, 42.43308162652588, -80.87301998715093],
# [126.08091298758552, 40.446149118689995, -82.61660874691425],
# [124.92210937602165, 43.59850481341899, -80.57809033527022],
# [123.20886721561682, 49.280519949260274, -81.88895195659163],
# [124.02163244359753, 50.102135196927755, -75.60177454432332],
# [124.67671288619638, 51.924743740595616, -78.52470196061927],
# [126.2969589612874, 52.207382607857774, -78.3407296179445],
# [127.51327332804047, 52.06341767016056, -77.5376573325906],
# [125.37353219920581, 52.45354184547641, -78.44345487465385],
# [123.5431571822743, 53.20902956285272, -77.89274730494311],
# [122.06613475735881, 53.656500065114784, -79.10815201933848],
# [121.22887581748051, 53.66822982429118, -80.15430849836198],
# [121.29863187132509, 53.057603616291345, -81.87794523343968],
# [121.47640281624213, 52.028962709642, -83.17853569530887],
# [120.36780816116726, 52.32157581157742, -84.29447086341364],
# [121.90107121222538, 53.150530854220996, -83.87516754041327],
# [122.30792162821903, 52.2225315441645, -82.82821783924972],
# [120.33177090353823, 57.690145508043166, -79.84526136454834],
# [120.38820843769486, 54.84868904002656, -79.2284333805261],
# [121.0662126535392, 52.099563698669634, -80.31673101244193],
# [122.79811327987035, 52.06689066593925, -79.80389361859554],
# [126.72341038019061, 48.44197536602091, -74.14489899045087],
# [126.10070552354854, 49.83499567705245, -76.56762024540362],
# [122.41888806373257, 51.43690412257744, -76.98792858319574],
# [120.24435440528023, 50.452765860108954, -78.37798597629923],
# [120.57281150541466, 52.05452098742229, -80.53833156368242],
# [121.90373678316102, 53.23271721388955, -80.50881309119073],
# [122.50255227291557, 53.27732244156752, -80.95864635228898],
# [121.87472416427146, 52.13787774272278, -82.22521992655209],
# [120.19951998500719, 49.957654228353654, -84.51289237325054],
# [120.20945526930386, 50.99029033194172, -85.57865463050086],
# [120.83621720839857, 50.75772970017323, -85.80497453133336],
# [121.924606218236, 53.29587274453488, -83.12484088210623],
# [122.74500723833663, 53.5329601704036, -82.24644736773037],
# [121.90952190731535, 53.98719313616031, -80.40074867955443],
# [122.74747256002446, 51.53239038745044, -80.75251377472372],
# [122.81482972194569, 51.515952259583436, -81.7373402966515],
# [122.82027972017897, 52.183989968038006, -82.72364151875102],
# [124.12104779800448, 51.17263510852614, -82.02362385914799],
# [123.85704137245078, 50.59759802057925, -80.83468732483263],
# [122.71095170989784, 47.93567401507012, -82.2220919937175],
# [123.02708199078486, 44.975076381512935, -84.34960435124269],
# [122.12730447736531, 45.501860082322594, -84.84324441877239],
# [122.78292996708528, 46.21554948246067, -84.67659132971438],
# [121.77881819184, 44.11585413772487, -84.06143150455806],
# [118.12198388061059, 45.11951056872483, -82.81445532934798],
# [113.17372560039337, 50.58043386484639, -86.51985540630852],
# [114.79621896475246, 50.86229086271987, -88.83448915792205],
# [121.06816949739446, 57.222230434366516, -79.94436036214718],
# [126.72735701258615, 56.39027631097382, -78.11876576948569],
# [130.07992237192263, 48.1605715700951, -78.16078004057042],
# [126.14350560119958, 43.84770791319735, -77.71689733464956],
# [123.70902848137058, 44.29434526219234, -78.94982052458053],
# [124.05812866406389, 47.442356735972666, -79.46219098402622],
# [123.01460734721469, 47.69268604494444, -81.88524671044541],
# [120.87922680395309, 50.49707623777617, -85.53920513061806],
# [120.3614087553321, 50.833837001720426, -84.95842885150174],
# [119.64029248140082, 50.9430142592573, -84.41478230610382],
# [121.78108262597391, 49.72443894204239, -80.43985053835443],
# [124.26270073667364, 46.44356426825039, -77.10363257374874],
# [124.18407493638902, 47.212840530544554, -77.32537358100535],
# [124.6156445146775, 47.866024187019455, -76.88524733122604],
# [124.35930018326616, 50.77998885806044, -79.11456581265016],
# [124.65673007419858, 47.710152591917414, -80.26546172013038],
# [123.96537037989192, 46.11889470260173, -81.16559744114622],
# [126.38793033072459, 46.7137623088924, -79.06334955408236],
# [129.2014924184205, 46.539301573671594, -75.01644888184451],
# [129.4284452599164, 47.57277911856106, -75.47456818478689],
# [128.42541321561137, 47.87244279888288, -76.25872484401768],
# [126.13014050388966, 47.53606824943342, -79.48739074836995],
# [126.34559043577332, 49.37005700743438, -80.91611244403887],
# [125.89553464179953, 49.75888293260284, -81.91811744773534],
# [124.42450986307261, 48.5994133191948, -84.77598790685316],
# [123.80930741710397, 46.72868050249324, -83.49132635838907],
# [123.92290681066554, 46.07689011284601, -82.97893541680683],
# [122.63356200196837, 47.3961654079915, -82.73538781336242],
# [124.11673832938345, 47.373024236874215, -82.87037879799195],
# [123.84659071180712, 49.456282066290086, -87.56123398888938],
# [123.1345448771306, 54.125597769102775, -88.25880462294631],
# [122.24336799304238, 54.24726640449067, -86.80616652702477],
# [127.66425124312673, 44.92255361121112, -72.86682224904274],
# [132.7600312113602, 38.52068652465863, -68.69162192412064],
# [136.7491229411954, 10.98558542674471, -73.42669146794265],
# [141.44277302010764, 4.862091910316131, -62.62563277966418],
# [132.37174676978594, 24.62604458433697, -85.77319287235233],
# [122.28751736074358, 48.74451335341505, -89.96990848393278],
# [123.44656743876115, 50.268951636142596, -85.8554677613021],
# [121.80395605698709, 57.98227912488444, -88.73274899331133],
# [121.24513731391008, 57.04363690811577, -89.91483037524957],
# [124.81735825857228, 50.128188865079835, -84.39364556668879],
# [125.15400045476618, 50.66661274122844, -77.77985830323753],
# [125.87008780249079, 44.04602714692647, -83.75895209516003],
# [126.32047995592778, 44.96270049990073, -82.06759385171456],
# [124.29238867739927, 43.015774638354785, -83.14462251808709],
# [124.8631052425201, 44.58165016216704, -83.90457298260088],
# [123.66742485951838, 46.83779961385656, -84.24895851113278],
# [119.54251530317403, 49.38665167965031, -84.51549616357646],
# [122.37082483097487, 53.1547788966909, -81.19598556270422],
# [121.55991476772124, 55.37981162217143, -81.57978040909626],
# [121.96957461309243, 55.09744613989357, -81.60372521731031],
# [121.93365895434397, 54.2784692333714, -83.017414776186],
# [121.55163160043791, 50.80958620685768, -80.9913336569055],
# [123.27704105006073, 46.88612757488302, -77.6852071313329],
# [126.86062348205587, 44.837932873769574, -79.23379704763425],
# [128.3598772953566, 47.480412391676374, -80.03588208736085],
# [127.34045346364555, 51.02509681837223, -77.68955192155133],
# [127.63389020644557, 51.639427309740405, -73.51399353574205],
# [129.29845805171627, 51.260438179038026, -72.76735452065702],
# [127.43680553204634, 48.57111977629498, -73.6724939662554],
# [129.25704063068937, 47.99519061378621, -75.43954358855912],
# [128.70835101376696, 48.00981260069866, -80.04129680460328],
# [128.5308895031276, 48.247269909415074, -80.54073981446415],
# [128.1605249688411, 47.65639103611723, -80.54593151896633],
# [128.28046759748642, 47.180925607773126, -79.90826572380263],
# [130.66270690888032, 44.717435250287025, -78.07737150066338],
# [129.8589427482224, 41.80123716091526, -78.38668850270761],
# [130.01487430740713, 42.4106400234136, -78.97274061299916],
# [126.69056027344861, 46.754376911753305, -78.9727657576152],
# [126.35362656068737, 49.232353232491675, -78.45947362841895],
# [127.82405213156034, 51.09304607331604, -76.078679358608],
# [128.90837307035315, 47.92321784113386, -76.50140324686059],
# [130.2847889125774, 48.05361413180153, -77.65351029568265],
# [128.39085553788837, 51.56519215107559, -77.45302178033864],
# [124.31505716304419, 51.92155634020901, -79.2068045840693],
# [118.30224775216625, 53.98670655723965, -79.65044859342706],
# [119.55091803850388, 55.55756074520639, -80.92354262634608],
# [118.45915069048154, 55.187495250946526, -84.57353050273146],
# [118.37386456360386, 53.527223529720466, -83.9340852491464],
# [119.5255226858814, 53.587449804475206, -82.15908445114232],
# [118.16766544196739, 51.97376369893894, -80.90787370939644],
# [117.30047054395313, 51.97191991436604, -80.60498459580137],
# [118.94774609261974, 51.9386129464815, -82.99001316333526],
# [120.45884858688983, 52.66369625102837, -85.64043675570005],
# [122.58540494363348, 50.370934595725075, -83.26534641876773],
# [125.66057133149727, 47.97382039777597, -81.40530101576076],
# [126.23225640553049, 47.786310354016095, -82.1339271968971],
# [126.36038384669398, 47.42051342380758, -84.04920114541167],
# [124.60691399342927, 49.46033930593012, -85.20379831291652],
# [120.7333639066095, 50.58614879865567, -82.32811371441493],
# [118.14218967741009, 52.80694997556333, -78.95470974280094],
# [119.22351973185297, 52.82610434139876, -79.32779980391312],
# [117.04895860638479, 50.418258294945566, -80.58747580294525],
# [118.79168091555607, 49.35092650937914, -81.31801680693769],
# [121.39660299790525, 47.64328166400309, -81.48609548355564],
# [122.13607296339747, 47.04568480049255, -80.22218926411549],
# [122.51470569701496, 48.31911312991832, -79.90277494612445],
# [121.66251057939344, 47.152350516809165, -80.78715753956435],
# [119.58971549124614, 48.578742652517874, -81.82789505139745],
# [118.01602132506345, 49.673553436604635, -83.3053699375656],
# [118.72946275332555, 51.55450474970651, -84.45813613604973],
# [120.82969284387184, 50.959127479489624, -85.43947683239051],
# [120.2206532259603, 51.728740584781576, -85.98312999266787],
# [121.51399714948411, 49.562484182039526, -84.10123395763289],
# [119.70978777189839, 49.46753006613002, -79.94318246912533],
# [119.07174187740543, 50.26247214895557, -78.82803378447362],
# [120.88499403354054, 51.66241160742615, -81.51344155520097],
# [121.87095146484918, 51.394581088255904, -83.61788670502452],
# [121.62339854236983, 47.2936470227532, -86.81482793078685],
# [121.00757587155194, 45.502461550136665, -92.41208585639048],
# [123.3343193213712, 64.94413525707019, -98.79182271523422],
# [-139.44444236806208, -70.133617697163, 87.83243879832463],
# [-143.7770730825155, -61.311285365741, 81.46484944917492],
# [-147.75808956394357, -62.99757478136102, 80.19838994012963],
# [138.92935709747093, 58.13799171147366, -83.31402838300986],
# [135.22979577378885, 53.394487737095375, -84.05016948458753],
# [136.5701725486534, 50.96314746322439, -77.69188076422837],
# [135.84485348613313, 53.89730019638694, -73.33694009302204],
# [136.74452050649418, 53.390351776259195, -70.56383930887672],
# [137.46185367115052, 53.13983044388964, -69.45242702501581],
# [137.48689685761045, 47.205212157040535, -73.46460376566061],
# [138.67093672153212, 46.52133946662025, -73.76355890551689],
# [142.52885233801237, 45.32965088361381, -70.70797824104595],
# [144.41753982413113, 46.47599105348469, -65.5364021802054],
# [147.98085934477, 50.54165659635774, -66.33627773993182],
# [151.8588931896769, 50.09751673399756, -66.84340838070891],
# [150.83096325664005, 45.89864913620907, -64.00584006681775],
# [147.1456677760226, 37.98166090590489, -66.68881799551791],
# [142.9788843424812, 35.980212408102254, -68.58687606044069],
# [143.83965049475273, 36.578809521006534, -69.20107432235832],
# [143.62848472926783, 36.61167628698464, -69.56219684620145],
# [143.65098878823397, 36.36381198887542, -70.19132612844732],
# [144.4273508863378, 36.82443320387479, -69.01711250186537],
# [144.37612284744014, 36.91114576044101, -69.07494630280571],
# [144.67785419282714, 37.0568355050116, -68.81934350009566],
# [147.66948178155434, 37.767663814197185, -64.75557010783653],
# [148.2110476508122, 34.39392191985721, -62.95918016102645],
# [152.10964286891112, 31.346970217727115, -61.628691583560624],
# [153.86837376884887, 31.976966489488603, -60.84747205846237],
# [154.52501262170736, 33.601236221754746, -67.46067098625475],
# [154.30236765739906, 32.245522607264284, -72.85110442932401],
# [158.54321153496494, 31.625265880081557, -70.2140075390069],
# [160.59302825171054, 34.41366500656412, -70.76017642162724],
# [-157.75001987907424, -36.697864003572214, 74.71431386868908],
# [-155.49691345554177, -37.015630876472926, 78.91510047050001],
# [-153.76715734850347, -39.337604809556936, 71.92679828160203],
# [-148.81576670807098, -44.5452154215307, 67.55934926967466],
# [-154.54712185552668, -48.45315314946549, 67.77798546345817],
# [-153.94770428218294, -59.915201731648814, 68.09298314309981],
# [-156.97455218578494, -55.78302052200982, 67.45349792395555],
# [153.4202668669928, 39.81427532166545, -59.261355920285865],
# [148.17129783285398, 28.71070011626827, -56.0461672700333],
# [148.04509602276053, 35.00522064229833, -56.541167623245386],
# [144.01099429517618, 38.97260841665631, -61.0213191848575],
# [137.7970610000713, 36.65847891201044, -68.77991698572384],
# [139.99642406290312, 33.19152999929201, -65.58176706052994],
# [146.66466995894513, 27.613458721685802, -61.093487965875234],
# [150.86361008353035, 23.479852280550435, -57.56410947118199],
# [153.58520271178028, 29.81243529563411, -55.47740076585592],
# [155.30783138106577, 32.06753275252949, -61.95066724910286],
# [159.12084380063945, 30.139930291123893, -66.52062316390136],
# [150.51461341669295, 32.746178797931044, -68.70544753945869],
# [142.17178133580393, 29.169567501466776, -75.65282257493683],
# [141.17870821154685, 30.83265290411937, -77.77886915139312],
# [141.9301199140175, 37.84698616845015, -74.82809656002465],
# [141.0381362383797, 39.259472521239665, -74.6399137672984],
# [138.75286178072017, 41.89482252441044, -72.48658484300823],
# [137.97480625522368, 44.196882473891, -71.62189767001861],
# [137.69329523960192, 45.85698379670635, -73.75680935896041],
# [133.7336624616985, 49.458798547140944, -76.83856606620142],
# [133.5586744416121, 49.57634840197156, -78.00908142942106],
# [132.37870590551984, 46.13060710932738, -78.39316755784982],
# [131.95836312747772, 43.03451507869912, -78.3214313334993],
# [132.3604051982243, 42.304233957133114, -78.88991795564633],
# [132.1915459462659, 43.42521275276891, -79.54611194310989],
# [132.9871627713596, 44.954729610694976, -80.21606276475526],
# [136.08457785378656, 45.01976943206674, -76.35629386406863],
# [137.66253868645035, 42.93326238456613, -74.46867758381006],
# [137.1399533808342, 46.43872449813158, -74.47799902572908],
# [135.40460927765034, 49.09529466538562, -74.52284110225084],
# [134.43497161571298, 48.09757159361239, -75.34815558006761],
# [134.6416880366954, 50.16457972750968, -77.4172781627294],
# [134.93763905794566, 48.84355848160631, -77.35701904961478],
# [136.97936300026714, 46.177117630279604, -77.55736922206324],
# [136.39905574333213, 42.78975803591761, -76.80567083689307],
# [137.75731766122024, 39.28104818479096, -76.81395268525021],
# [141.2642188394449, 37.60554091292037, -76.31164235178645],
# [143.01635251473155, 35.795184728619766, -73.548877955561],
# [143.2171228982183, 34.37511656004443, -73.6089870637451],
# [145.52885474951114, 31.663176764698488, -72.87318615782405],
# [147.4499221155343, 30.76006990502134, -68.94842977248878],
# [147.6685827859223, 34.85870151780462, -66.72632878385677],
# [147.7352618143689, 41.25334490344059, -67.03077240113643],
# [146.65751354932823, 42.672735440859505, -68.09580138274268],
# [147.49485620510427, 43.93434784791752, -67.19988745129267],
# [148.11986221257732, 48.80107968686981, -65.35031065289503],
# [149.86214983361953, 50.58903905493068, -69.22554229411843],
# [155.23430177075602, 48.466105869530175, -69.04437050607196],
# [156.20564432307555, 51.2746735302934, -66.71744779640514],
# [-152.1528189822149, -64.04411364913835, 66.3288895547715],
# [141.34254149605027, 65.66883850537953, -76.50004067466615],
# [116.4287646101505, 79.0057725004109, -62.152453989635184],
# [109.95347191512144, 82.65018731060484, -76.02976035812102],
# [128.63942060619527, 35.72099671951293, -67.98080219408945],
# [132.6725338096367, 19.15524641886345, -59.120666486929714],
# [145.20348351830486, 17.91032458445531, -67.3365727962],
# [144.73630641924268, 21.83335104584427, -70.94111118970679],
# [140.76639634729688, 29.668150517580543, -77.86770766953676],
# [141.891052407311, 41.20479146762281, -75.68881123772414],
# [138.7917094260748, 41.1845174043082, -77.35092501142519],
# [135.70270533551056, 42.553495850336866, -75.38569663318344],
# [131.32677565059558, 43.7821246321273, -74.53752746383049],
# [132.16942859682274, 43.25457353762908, -73.97017959612934],
# [135.03786001872575, 40.470692102977054, -75.55588303321457],
# [133.2879235669383, 38.132519421477035, -79.24881477321544],
# [133.49356727728434, 37.230002106759294, -80.15769905188746],
# [131.46761707725796, 38.49303118437562, -79.79076293534459],
# [130.5269020497727, 41.722823174869106, -76.36184074645581],
# [130.04207328590016, 41.78300264056538, -77.70361101930244],
# [129.12099348034255, 42.15916137323013, -76.71130640708104],
# [126.8251064373025, 42.281002272843864, -75.5834449699478],
# [126.54961510968668, 44.255581208954624, -74.51905581774349],
# [124.91439861413735, 45.15407995182321, -74.46035837148675],
# [125.04324714626651, 48.921583588119276, -79.03970903544476],
# [123.98327451614985, 50.669267364430254, -79.73107637634105],
# [122.01085914680215, 47.17927774811161, -78.74536741498721],
# [126.67912383321419, 47.62930678469444, -80.63889790433633],
# [126.40022062242912, 48.26758472816983, -79.2697991565744],
# [126.84025980695299, 49.792993329300344, -80.45699662560843],
# [123.54947249766047, 52.51026201159924, -83.93271375606676],
# [119.71965916267236, 53.32197170020214, -83.78255361742423],
# [120.2477222331202, 51.299063869397166, -83.51274554539103],
# [119.83846458735778, 49.50939560221172, -81.52064978170546],
# [122.14240196302731, 48.51844921192265, -82.29190868659288],
# [121.82792359531095, 49.0380025811688, -81.80625621440936],
# [121.28918881512043, 47.587342758424626, -80.73315567340732],
# [120.42654550052121, 48.498694455115036, -83.32060971020329],
# [124.3236327163625, 48.42162246819251, -83.09153978369986],
# [126.85655682269646, 48.89730443991738, -83.70249606031294],
# [127.83934912008849, 50.43581672239281, -82.20009235126075],
# [127.55530234042375, 50.54742592433999, -78.90607929032814],
# [125.68168064353937, 48.49464346307414, -79.0659434870821],
# [126.52065731271281, 46.13114241745373, -80.59990674931603],
# [126.32154087926534, 47.17739510137365, -81.17769416337376],
# [126.15861413330728, 47.48767311358681, -80.49135398042785],
# [124.64359256259816, 50.74311105238352, -77.52814398040942],
# [123.91155210192328, 51.953348080998744, -78.91815327310084],
# [121.22088216553395, 52.28533935927183, -81.04721393266537],
# [118.42177269154067, 50.43823916131004, -81.84789741238619],
# [119.29775106266617, 47.54110507729978, -84.26479419195722],
# [118.84599638181707, 48.40434767630053, -84.26026342398637],
# [118.90548094777435, 48.51239048334802, -83.47529013403522],
# [121.96616201790941, 50.5299721380281, -84.52378018704457],
# [121.73893033142181, 52.46427653512619, -86.7046353750173],
# [121.16650826569852, 47.789899528647766, -84.65579673567414],
# [122.18601018848967, 45.45961138188503, -82.48915362374136],
# [121.69449277939552, 44.637474044928254, -82.2719074928137],
# [120.98443760679568, 46.06035746524875, -83.25497301815217],
# [120.44424070910038, 47.41582857035238, -81.09421132076133],
# [118.97571340624607, 50.158352619112186, -81.86300620435968],
# [121.0788228288074, 57.9495490520442, -85.76427021983825],
# [120.88714656149173, 52.59898742688067, -85.80499484916848],
# [120.60389106292187, 46.890014443558854, -87.02596724948053],
# [118.86882914042313, 46.7207879903393, -86.21867401754822],
# [113.01203693995943, 48.05672439349338, -87.91291985592859],
# [116.336379160459, 46.24640731312869, -89.09708881093272],
# [115.56169610859023, 44.68951820832176, -93.18853487591122],
# [114.26919211625527, 41.426145738115615, -94.59368683441056],
# [116.22493418680338, 39.174212921932025, -93.47361878135105],
# [117.55052160576004, 38.09128701342917, -93.48261190002857],
# [118.68629708989208, 38.40043946613556, -92.51133543199494],
# [120.27595528427686, 38.67019206376213, -92.51675629899273]])
#
# # Define the number of clusters (k)
# k = 3
#
# # Initialize cluster centroids randomly
# centroids = X[np.random.choice(X.shape[0], k, replace=False)]
#
# # Maximum number of iterations
# max_iters = 100
#
# # K-means algorithm
# for _ in range(max_iters):
#     # Assign each data point to the nearest cluster
#     distances = np.linalg.norm(X[:, np.newaxis, :] - centroids, axis=2)
#     labels = np.argmin(distances, axis=1)
#
#     # Update cluster centroids
#     new_centroids = np.array([X[labels == i].mean(axis=0) for i in range(k)])
#
#     # Check for convergence
#     if np.all(centroids == new_centroids):
#         break
#
#     centroids = new_centroids
#
# # Final cluster assignments
# print(labels)
# unique, counts = np.unique(labels, return_counts=True)
# print(f"Unique = {unique}")
# print(f"Counts = {counts}")

