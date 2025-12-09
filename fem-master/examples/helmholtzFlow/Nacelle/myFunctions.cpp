#include "gmsh.h"

#include <fstream>
#include <gmshddm/Formulation.h>
#include <gmshddm/GmshDdm.h>
#include <gmshfem/AnalyticalFunction.h>
#include <gmshfem/Message.h>
#include <iostream>

using namespace gmshddm;
using namespace gmshddm::common;
using namespace gmshddm::domain;
using namespace gmshddm::problem;
using namespace gmshddm::field;

using namespace gmshfem;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::analytics;
using namespace gmshfem::post;
using namespace gmshfem::equation;

// ***************************************************************************************

void Nacelle(const double h_near, const double h_far, const double Lx, const double Ly, const double L_duct, int ElementOrder, int Npml, int Npml_active)
{
  gmsh::option::setNumber("General.Terminal", 1);
  gmsh::model::add("Nacelle");
  gmsh::option::setNumber("Mesh.CharacteristicLengthMin", h_near);
  gmsh::option::setNumber("Mesh.CharacteristicLengthMax", h_far);

  // Nacelle (points)
  gmsh::model::occ::addPoint(L_duct, 1.5000000000000000, 0, h_near, 1);
  gmsh::model::occ::addPoint(-0.2987304221999528, 1.4999980867903138, 0, h_near, 2);
  gmsh::model::occ::addPoint(-0.2954590219195146, 1.4999894098722784, 0, h_near, 3);
  gmsh::model::occ::addPoint(-0.2902590561329333, 1.4999695633124288, 0, h_near, 4);
  gmsh::model::occ::addPoint(-0.2832037818144567, 1.4999341411772997, 0, h_near, 5);
  gmsh::model::occ::addPoint(-0.2743664559383324, 1.4998787375334264, 0, h_near, 6);
  gmsh::model::occ::addPoint(-0.2638203354788085, 1.4997989464473434, 0, h_near, 7);
  gmsh::model::occ::addPoint(-0.2516386774101325, 1.4996903619855857, 0, h_near, 8);
  gmsh::model::occ::addPoint(-0.2378947387065523, 1.4995485782146885, 0, h_near, 9);
  gmsh::model::occ::addPoint(-0.2226617763423157, 1.4993691892011862, 0, h_near, 10);
  gmsh::model::occ::addPoint(-0.2060130472916705, 1.4991477890116141, 0, h_near, 11);
  gmsh::model::occ::addPoint(-0.1880218085288645, 1.4988799717125068, 0, h_near, 12);
  gmsh::model::occ::addPoint(-0.1687613170281455, 1.4985613313703996, 0, h_near, 13);
  gmsh::model::occ::addPoint(-0.1483048297637612, 1.4981874620518270, 0, h_near, 14);
  gmsh::model::occ::addPoint(-0.1267256037099594, 1.4977539578233241, 0, h_near, 15);
  gmsh::model::occ::addPoint(-0.1040968958409880, 1.4972564127514256, 0, h_near, 16);
  gmsh::model::occ::addPoint(-0.0804919631310948, 1.4966904209026668, 0, h_near, 17);
  gmsh::model::occ::addPoint(-0.0559840625545274, 1.4960515763435822, 0, h_near, 18);
  gmsh::model::occ::addPoint(-0.0306464510855338, 1.4953354731407069, 0, h_near, 19);
  gmsh::model::occ::addPoint(-0.0045523856983616, 1.4945377053605760, 0, h_near, 20);
  gmsh::model::occ::addPoint(0.0222248766327412, 1.4936543132180056, 0, h_near, 21);
  gmsh::model::occ::addPoint(0.0496120789335270, 1.4926865234015851, 0, h_near, 22);
  gmsh::model::occ::addPoint(0.0775359642297478, 1.4916389087120159, 0, h_near, 23);
  gmsh::model::occ::addPoint(0.1059232755471561, 1.4905160977185339, 0, h_near, 24);
  gmsh::model::occ::addPoint(0.1347007559115039, 1.4893227189903755, 0, h_near, 25);
  gmsh::model::occ::addPoint(0.1637951483485434, 1.4880634010967773, 0, h_near, 26);
  gmsh::model::occ::addPoint(0.1931331958840269, 1.4867427726069748, 0, h_near, 27);
  gmsh::model::occ::addPoint(0.2226416415437066, 1.4853654620902048, 0, h_near, 28);
  gmsh::model::occ::addPoint(0.2522472283533346, 1.4839360981157035, 0, h_near, 29);
  gmsh::model::occ::addPoint(0.2818766993386633, 1.4824593092527072, 0, h_near, 30);
  gmsh::model::occ::addPoint(0.3114567975254448, 1.4809397240704518, 0, h_near, 31);
  gmsh::model::occ::addPoint(0.3409142659394311, 1.4793819711381737, 0, h_near, 32);
  gmsh::model::occ::addPoint(0.3701758476063748, 1.4777906790251094, 0, h_near, 33);
  gmsh::model::occ::addPoint(0.3991682855520280, 1.4761704763004948, 0, h_near, 34);
  gmsh::model::occ::addPoint(0.4278183228021426, 1.4745259915335662, 0, h_near, 35);
  gmsh::model::occ::addPoint(0.4560527023824713, 1.4728618532935602, 0, h_near, 36);
  gmsh::model::occ::addPoint(0.4837981673187660, 1.4711826901497127, 0, h_near, 37);
  gmsh::model::occ::addPoint(0.5109814606367789, 1.4694931306712600, 0, h_near, 38);
  gmsh::model::occ::addPoint(0.5375293253622622, 1.4677978034274384, 0, h_near, 39);
  gmsh::model::occ::addPoint(0.5633692282590527, 1.4661012625652645, 0, h_near, 40);
  gmsh::model::occ::addPoint(0.5884720603760486, 1.4644035968985738, 0, h_near, 41);
  gmsh::model::occ::addPoint(0.6128760204039941, 1.4626979739747705, 0, h_near, 42);
  gmsh::model::occ::addPoint(0.6366250969383089, 1.4609769659635023, 0, h_near, 43);
  gmsh::model::occ::addPoint(0.6597632785744117, 1.4592331450344158, 0, h_near, 44);
  gmsh::model::occ::addPoint(0.6823345539077219, 1.4574590833571586, 0, h_near, 45);
  gmsh::model::occ::addPoint(0.7043829115336588, 1.4556473531013774, 0, h_near, 46);
  gmsh::model::occ::addPoint(0.7259523400476415, 1.4537905264367195, 0, h_near, 47);
  gmsh::model::occ::addPoint(0.7470868280450892, 1.4518811755328322, 0, h_near, 48);
  gmsh::model::occ::addPoint(0.7678303641214211, 1.4499118725593625, 0, h_near, 49);
  gmsh::model::occ::addPoint(0.7882269368720566, 1.4478751896859574, 0, h_near, 50);
  gmsh::model::occ::addPoint(0.8083205348924145, 1.4457636990822642, 0, h_near, 51);
  gmsh::model::occ::addPoint(0.8281551467779145, 1.4435699729179299, 0, h_near, 52);
  gmsh::model::occ::addPoint(0.8477747611239757, 1.4412865833626016, 0, h_near, 53);
  gmsh::model::occ::addPoint(0.8672233665260169, 1.4389061025859267, 0, h_near, 54);
  gmsh::model::occ::addPoint(0.8865449515794576, 1.4364211027575522, 0, h_near, 55);
  gmsh::model::occ::addPoint(0.9057835048797171, 1.4338241560471250, 0, h_near, 56);
  gmsh::model::occ::addPoint(0.9249830150222146, 1.4311078346242925, 0, h_near, 57);
  gmsh::model::occ::addPoint(0.9441874706023690, 1.4282647106587019, 0, h_near, 58);
  gmsh::model::occ::addPoint(0.9634408602155999, 1.4252873563199999, 0, h_near, 59);
  gmsh::model::occ::addPoint(0.9827722478565070, 1.4221712172508760, 0, h_near, 60);
  gmsh::model::occ::addPoint(1.0021509991164135, 1.4189232329861872, 0, h_near, 61);
  gmsh::model::occ::addPoint(1.0215315549858226, 1.4155532165338327, 0, h_near, 62);
  gmsh::model::occ::addPoint(1.0408683564552379, 1.4120709809017116, 0, h_near, 63);
  gmsh::model::occ::addPoint(1.0601158445151637, 1.4084863390977229, 0, h_near, 64);
  gmsh::model::occ::addPoint(1.0792284601561029, 1.4048091041297661, 0, h_near, 65);
  gmsh::model::occ::addPoint(1.0981606443685594, 1.4010490890057399, 0, h_near, 66);
  gmsh::model::occ::addPoint(1.1168668381430371, 1.3972161067335436, 0, h_near, 67);
  gmsh::model::occ::addPoint(1.1353014824700389, 1.3933199703210764, 0, h_near, 68);
  gmsh::model::occ::addPoint(1.1534190183400692, 1.3893704927762371, 0, h_near, 69);
  gmsh::model::occ::addPoint(1.1711738867436312, 1.3853774871069253, 0, h_near, 70);
  gmsh::model::occ::addPoint(1.1885205286712286, 1.3813507663210400, 0, h_near, 71);
  gmsh::model::occ::addPoint(1.2054133851133650, 1.3773001434264802, 0, h_near, 72);
  gmsh::model::occ::addPoint(1.2218068970605445, 1.3732354314311450, 0, h_near, 73);
  gmsh::model::occ::addPoint(1.2376555055032699, 1.3691664433429336, 0, h_near, 74);
  gmsh::model::occ::addPoint(1.2529136514320454, 1.3651029921697451, 0, h_near, 75);
  gmsh::model::occ::addPoint(1.2675357758373744, 1.3610548909194788, 0, h_near, 76);
  gmsh::model::occ::addPoint(1.2814763197097605, 1.3570319526000336, 0, h_near, 77);
  gmsh::model::occ::addPoint(1.2946897240397075, 1.3530439902193088, 0, h_near, 78);
  gmsh::model::occ::addPoint(1.3071340397340678, 1.3491002687496183, 0, h_near, 79);
  gmsh::model::occ::addPoint(1.3188092829772511, 1.3452036822495976, 0, h_near, 80);
  gmsh::model::occ::addPoint(1.3297425443262842, 1.3413530145109940, 0, h_near, 81);
  gmsh::model::occ::addPoint(1.3399613655777376, 1.3375469808211058, 0, h_near, 82);
  gmsh::model::occ::addPoint(1.3494932885281816, 1.3337842964672317, 0, h_near, 83);
  gmsh::model::occ::addPoint(1.3583658549741870, 1.3300636767366700, 0, h_near, 84);
  gmsh::model::occ::addPoint(1.3666066067123239, 1.3263838369167200, 0, h_near, 85);
  gmsh::model::occ::addPoint(1.3742430855391632, 1.3227434922946795, 0, h_near, 86);
  gmsh::model::occ::addPoint(1.3813028332512749, 1.3191413581578477, 0, h_near, 87);
  gmsh::model::occ::addPoint(1.3878133916452300, 1.3155761497935230, 0, h_near, 88);
  gmsh::model::occ::addPoint(1.3938023025175985, 1.3120465824890037, 0, h_near, 89);
  gmsh::model::occ::addPoint(1.3992971076649514, 1.3085513715315891, 0, h_near, 90);
  gmsh::model::occ::addPoint(1.4043253488838587, 1.3050892322085772, 0, h_near, 91);
  gmsh::model::occ::addPoint(1.4089145679708910, 1.3016588798072668, 0, h_near, 92);
  gmsh::model::occ::addPoint(1.4130923067226191, 1.2982590296149565, 0, h_near, 93);
  gmsh::model::occ::addPoint(1.4168861069356131, 1.2948883969189451, 0, h_near, 94);
  gmsh::model::occ::addPoint(1.4203235104064436, 1.2915456970065309, 0, h_near, 95);
  gmsh::model::occ::addPoint(1.4234320589316809, 1.2882296451650126, 0, h_near, 96);
  gmsh::model::occ::addPoint(1.4262392943078961, 1.2849389566816889, 0, h_near, 97);
  gmsh::model::occ::addPoint(1.4287726303873629, 1.2816723494857822, 0, h_near, 98);
  gmsh::model::occ::addPoint(1.4310518043646006, 1.2784287000219530, 0, h_near, 99);
  gmsh::model::occ::addPoint(1.4330846546146057, 1.2752071304337900, 0, h_near, 100);
  gmsh::model::occ::addPoint(1.4348779959580076, 1.2720067840002738, 0, h_near, 101);
  gmsh::model::occ::addPoint(1.4364386432154359, 1.2688268040003843, 0, h_near, 102);
  gmsh::model::occ::addPoint(1.4377734112075198, 1.2656663337131024, 0, h_near, 103);
  gmsh::model::occ::addPoint(1.4388891147548886, 1.2625245164174086, 0, h_near, 104);
  gmsh::model::occ::addPoint(1.4397925686781721, 1.2594004953922833, 0, h_near, 105);
  gmsh::model::occ::addPoint(1.4404905877979992, 1.2562934139167068, 0, h_near, 106);
  gmsh::model::occ::addPoint(1.4409899869349998, 1.2532024152696595, 0, h_near, 107);
  gmsh::model::occ::addPoint(1.4412975809098028, 1.2501266427301221, 0, h_near, 108);
  gmsh::model::occ::addPoint(1.4414201845430381, 1.2470652395770749, 0, h_near, 109);
  gmsh::model::occ::addPoint(1.4413646126553346, 1.2440173490894983, 0, h_near, 110);
  gmsh::model::occ::addPoint(1.4411376800673219, 1.2409821145463731, 0, h_near, 111);
  gmsh::model::occ::addPoint(1.4407462015996295, 1.2379586792266792, 0, h_near, 112);
  gmsh::model::occ::addPoint(1.4401969920728868, 1.2349461864093974, 0, h_near, 113);
  gmsh::model::occ::addPoint(1.4394968663077232, 1.2319437793735082, 0, h_near, 114);
  gmsh::model::occ::addPoint(1.4386526391247678, 1.2289506013979918, 0, h_near, 115);
  gmsh::model::occ::addPoint(1.4376711253446504, 1.2259657957618291, 0, h_near, 116);
  gmsh::model::occ::addPoint(1.4365591397880000, 1.2229885057440000, 0, h_near, 117);
  gmsh::model::occ::addPoint(1.4353209654502130, 1.2200181664371075, 0, h_near, 118);
  gmsh::model::occ::addPoint(1.4339507580257536, 1.2170553801882427, 0, h_near, 119);
  gmsh::model::occ::addPoint(1.4324401413838517, 1.2141010411581195, 0, h_near, 120);
  gmsh::model::occ::addPoint(1.4307807393937382, 1.2111560435074511, 0, h_near, 121);
  gmsh::model::occ::addPoint(1.4289641759246439, 1.2082212813969513, 0, h_near, 122);
  gmsh::model::occ::addPoint(1.4269820748457991, 1.2052976489873333, 0, h_near, 123);
  gmsh::model::occ::addPoint(1.4248260600264344, 1.2023860404393112, 0, h_near, 124);
  gmsh::model::occ::addPoint(1.4224877553357806, 1.1994873499135981, 0, h_near, 125);
  gmsh::model::occ::addPoint(1.4199587846430684, 1.1966024715709076, 0, h_near, 126);
  gmsh::model::occ::addPoint(1.4172307718175279, 1.1937322995719533, 0, h_near, 127);
  gmsh::model::occ::addPoint(1.4142953407283902, 1.1908777280774490, 0, h_near, 128);
  gmsh::model::occ::addPoint(1.4111441152448858, 1.1880396512481077, 0, h_near, 129);
  gmsh::model::occ::addPoint(1.4077687192362451, 1.1852189632446435, 0, h_near, 130);
  gmsh::model::occ::addPoint(1.4041607765716990, 1.1824165582277697, 0, h_near, 131);
  gmsh::model::occ::addPoint(1.4003119111204778, 1.1796333303581996, 0, h_near, 132);
  gmsh::model::occ::addPoint(1.3962137467518123, 1.1768701737966474, 0, h_near, 133);
  gmsh::model::occ::addPoint(1.3918579073349331, 1.1741279827038260, 0, h_near, 134);
  gmsh::model::occ::addPoint(1.3872360167390707, 1.1714076512404492, 0, h_near, 135);
  gmsh::model::occ::addPoint(1.3823396988334558, 1.1687100735672304, 0, h_near, 136);
  gmsh::model::occ::addPoint(1.3771617800741638, 1.1660363424023938, 0, h_near, 137);
  gmsh::model::occ::addPoint(1.3717090669893410, 1.1633898586952205, 0, h_near, 138);
  gmsh::model::occ::addPoint(1.3659973855084711, 1.1607755125763191, 0, h_near, 139);
  gmsh::model::occ::addPoint(1.3600427118843919, 1.1581982189959870, 0, h_near, 140);
  gmsh::model::occ::addPoint(1.3538610223699419, 1.1556628929045218, 0, h_near, 141);
  gmsh::model::occ::addPoint(1.3474682932179600, 1.1531744492522209, 0, h_near, 142);
  gmsh::model::occ::addPoint(1.3408805006812845, 1.1507378029893816, 0, h_near, 143);
  gmsh::model::occ::addPoint(1.3341136210127535, 1.1483578690663012, 0, h_near, 144);
  gmsh::model::occ::addPoint(1.3271836304652063, 1.1460395624332775, 0, h_near, 145);
  gmsh::model::occ::addPoint(1.3201065052914809, 1.1437877980406079, 0, h_near, 146);
  gmsh::model::occ::addPoint(1.3128982217444156, 1.1416074908385898, 0, h_near, 147);
  gmsh::model::occ::addPoint(1.3055747560768494, 1.1395035557775204, 0, h_near, 148);
  gmsh::model::occ::addPoint(1.2981520845416201, 1.1374809078076973, 0, h_near, 149);
  gmsh::model::occ::addPoint(1.2906461833915670, 1.1355444618794182, 0, h_near, 150);
  gmsh::model::occ::addPoint(1.2830730288795280, 1.1336991329429800, 0, h_near, 151);
  gmsh::model::occ::addPoint(1.2754485972583420, 1.1319498359486806, 0, h_near, 152);
  gmsh::model::occ::addPoint(1.2677888647808473, 1.1303014858468172, 0, h_near, 153);
  gmsh::model::occ::addPoint(1.2601098076998822, 1.1287589975876873, 0, h_near, 154);
  gmsh::model::occ::addPoint(1.2524274022682853, 1.1273272861215882, 0, h_near, 155);
  gmsh::model::occ::addPoint(1.2447574699229378, 1.1260112034409910, 0, h_near, 156);
  gmsh::model::occ::addPoint(1.2371065431432720, 1.1248118240687639, 0, h_near, 157);
  gmsh::model::occ::addPoint(1.2294667565246746, 1.1237243674498927, 0, h_near, 158);
  gmsh::model::occ::addPoint(1.2218290061348731, 1.1227435493667490, 0, h_near, 159);
  gmsh::model::occ::addPoint(1.2141841880415942, 1.1218640856017048, 0, h_near, 160);
  gmsh::model::occ::addPoint(1.2065231983125653, 1.1210806919371317, 0, h_near, 161);
  gmsh::model::occ::addPoint(1.1988369330155138, 1.1203880841554017, 0, h_near, 162);
  gmsh::model::occ::addPoint(1.1911162882181663, 1.1197809780388865, 0, h_near, 163);
  gmsh::model::occ::addPoint(1.1833521599882502, 1.1192540893699578, 0, h_near, 164);
  gmsh::model::occ::addPoint(1.1760000000000002, 1.1188021339309877, 0, h_near, 165);
  gmsh::model::occ::addPoint(1.1676570375016206, 1.1184198275043478, 0, h_near, 166);
  gmsh::model::occ::addPoint(1.1597078353803612, 1.1181018858724101, 0, h_near, 167);
  gmsh::model::occ::addPoint(1.1516787340974417, 1.1178430248175462, 0, h_near, 168);
  gmsh::model::occ::addPoint(1.1435606297205891, 1.1176379601221278, 0, h_near, 169);
  gmsh::model::occ::addPoint(1.1353444183175307, 1.1174814075685273, 0, h_near, 170);
  gmsh::model::occ::addPoint(1.1270209959559938, 1.1173680829391159, 0, h_near, 171);
  gmsh::model::occ::addPoint(1.1185812587037050, 1.1172927020162657, 0, h_near, 172);
  gmsh::model::occ::addPoint(1.1100161026283915, 1.1172499805823484, 0, h_near, 173);
  gmsh::model::occ::addPoint(1.1013164237977808, 1.1172346344197359, 0, h_near, 174);
  gmsh::model::occ::addPoint(1.0924731182795999, 1.1172413793108000, 0, h_near, 175);
  gmsh::model::occ::addPoint(1.0834748637432567, 1.1172672429155304, 0, h_near, 176);
  gmsh::model::occ::addPoint(1.0743014642648825, 1.1173185004043873, 0, h_near, 177);
  gmsh::model::occ::addPoint(1.0649305055222895, 1.1174037388254490, 0, h_near, 178);
  gmsh::model::occ::addPoint(1.0553395731932902, 1.1175315452267940, 0, h_near, 179);
  gmsh::model::occ::addPoint(1.0455062529556964, 1.1177105066565003, 0, h_near, 180);
  gmsh::model::occ::addPoint(1.0354081304873206, 1.1179492101626463, 0, h_near, 181);
  gmsh::model::occ::addPoint(1.0250227914659749, 1.1182562427933103, 0, h_near, 182);
  gmsh::model::occ::addPoint(1.0143278215694715, 1.1186401915965702, 0, h_near, 183);
  gmsh::model::occ::addPoint(1.0033008064756224, 1.1191096436205046, 0, h_near, 184);
  gmsh::model::occ::addPoint(0.9919193318622402, 1.1196731859131916, 0, h_near, 185);
  gmsh::model::occ::addPoint(0.9801609834071370, 1.1203394055227096, 0, h_near, 186);
  gmsh::model::occ::addPoint(0.9680033467881249, 1.1211168894971368, 0, h_near, 187);
  gmsh::model::occ::addPoint(0.9554240076830162, 1.1220142248845513, 0, h_near, 188);
  gmsh::model::occ::addPoint(0.9424005517696229, 1.1230399987330313, 0, h_near, 189);
  gmsh::model::occ::addPoint(0.9289105647257571, 1.1242027980906553, 0, h_near, 190);
  gmsh::model::occ::addPoint(0.9149316322292316, 1.1255112100055016, 0, h_near, 191);
  gmsh::model::occ::addPoint(0.9004413399578581, 1.1269738215256480, 0, h_near, 192);
  gmsh::model::occ::addPoint(0.8854172735894490, 1.1285992196991732, 0, h_near, 193);
  gmsh::model::occ::addPoint(0.8698370188018166, 1.1303959915741550, 0, h_near, 194);
  gmsh::model::occ::addPoint(0.8536802654824547, 1.1323714304349906, 0, h_near, 195);
  gmsh::model::occ::addPoint(0.8369511649564076, 1.1345177895632750, 0, h_near, 196);
  gmsh::model::occ::addPoint(0.8196696501213343, 1.1368176190129913, 0, h_near, 197);
  gmsh::model::occ::addPoint(0.8018559169011027, 1.1392533071176614, 0, h_near, 198);
  gmsh::model::occ::addPoint(0.7835301612195810, 1.1418072422108076, 0, h_near, 199);
  gmsh::model::occ::addPoint(0.7647125790006376, 1.1444618126259525, 0, h_near, 200);
  gmsh::model::occ::addPoint(0.7454233661681409, 1.1471994066966180, 0, h_near, 201);
  gmsh::model::occ::addPoint(0.7256827186459585, 1.1500024127563266, 0, h_near, 202);
  gmsh::model::occ::addPoint(0.7055108323579597, 1.1528532191386007, 0, h_near, 203);
  gmsh::model::occ::addPoint(0.6849279032280126, 1.1557342141769620, 0, h_near, 204);
  gmsh::model::occ::addPoint(0.6639541271799853, 1.1586277862049335, 0, h_near, 205);
  gmsh::model::occ::addPoint(0.6426097001377460, 1.1615163235560368, 0, h_near, 206);
  gmsh::model::occ::addPoint(0.6209148180251627, 1.1643822145637948, 0, h_near, 207);
  gmsh::model::occ::addPoint(0.5988896767661047, 1.1672078475617294, 0, h_near, 208);
  gmsh::model::occ::addPoint(0.5765544722844398, 1.1699756108833628, 0, h_near, 209);
  gmsh::model::occ::addPoint(0.5539294005040362, 1.1726678928622174, 0, h_near, 210);
  gmsh::model::occ::addPoint(0.5310346573487623, 1.1752670818318154, 0, h_near, 211);
  gmsh::model::occ::addPoint(0.5078904387424863, 1.1777555661256793, 0, h_near, 212);
  gmsh::model::occ::addPoint(0.4845169406090762, 1.1801157340773314, 0, h_near, 213);
  gmsh::model::occ::addPoint(0.4609343588724011, 1.1823301111174744, 0, h_near, 214);
  gmsh::model::occ::addPoint(0.4371628894563291, 1.1843894485076563, 0, h_near, 215);
  gmsh::model::occ::addPoint(0.4132227282847282, 1.1862972475472342, 0, h_near, 216);
  gmsh::model::occ::addPoint(0.3891340712814669, 1.1880581063130116, 0, h_near, 217);
  gmsh::model::occ::addPoint(0.3649171143704134, 1.1896766228817914, 0, h_near, 218);
  gmsh::model::occ::addPoint(0.3405920534754354, 1.1911573953303773, 0, h_near, 219);
  gmsh::model::occ::addPoint(0.3161790845204026, 1.1925050217355722, 0, h_near, 220);
  gmsh::model::occ::addPoint(0.2916984034291824, 1.1937241001741792, 0, h_near, 221);
  gmsh::model::occ::addPoint(0.2671702061256434, 1.1948192287230019, 0, h_near, 222);
  gmsh::model::occ::addPoint(0.2426146885336538, 1.1957950054588433, 0, h_near, 223);
  gmsh::model::occ::addPoint(0.2180520465770818, 1.1966560284585066, 0, h_near, 224);
  gmsh::model::occ::addPoint(0.1935024761797952, 1.1974068957987953, 0, h_near, 225);
  gmsh::model::occ::addPoint(0.1689861732656636, 1.1980522055565124, 0, h_near, 226);
  gmsh::model::occ::addPoint(0.1445233337585545, 1.1985965558084610, 0, h_near, 227);
  gmsh::model::occ::addPoint(0.1201341535823364, 1.1990445446314448, 0, h_near, 228);
  gmsh::model::occ::addPoint(0.0960000000000000, 1.1994007701022666, 0, h_near, 229);
  gmsh::model::occ::addPoint(0.0716575549180462, 1.1996698302977298, 0, h_near, 230);
  gmsh::model::occ::addPoint(0.0476105282777101, 1.1998563232946375, 0, h_near, 231);
  gmsh::model::occ::addPoint(0.0237179446637388, 1.1999648471697932, 0, h_near, 232);
  gmsh::model::occ::addPoint(0.0000000000000000, 1.2000000000000000, 0, h_near, 233);

  // Fan duct (points)
  gmsh::model::occ::addPoint(L_duct, 1.2000000000000000, 0, h_near, 234);
  gmsh::model::occ::addPoint(L_duct, 0.3586206896556000, 0, h_near, 235);

  // Spinner (points)
  gmsh::model::occ::addPoint(0.0000000000000000, 0.3586206896556000, 0, h_near, 236);
  gmsh::model::occ::addPoint(0.0094369356474370, 0.3585833677801259, 0, h_near, 237);
  gmsh::model::occ::addPoint(0.0194167690446929, 0.3584636577634160, 0, h_near, 238);
  gmsh::model::occ::addPoint(0.0299032528379802, 0.3582499430200392, 0, h_near, 239);
  gmsh::model::occ::addPoint(0.0408601396735112, 0.3579306069645641, 0, h_near, 240);
  gmsh::model::occ::addPoint(0.0522511821974983, 0.3574940330115595, 0, h_near, 241);
  gmsh::model::occ::addPoint(0.0640401330561540, 0.3569286045755942, 0, h_near, 242);
  gmsh::model::occ::addPoint(0.0761907448956906, 0.3562227050712369, 0, h_near, 243);
  gmsh::model::occ::addPoint(0.0886667703623206, 0.3553647179130562, 0, h_near, 244);
  gmsh::model::occ::addPoint(0.1014319621022565, 0.3543430265156212, 0, h_near, 245);
  gmsh::model::occ::addPoint(0.1144500727617104, 0.3531460192702524, 0, h_near, 246);
  gmsh::model::occ::addPoint(0.1276848549868951, 0.3517631396397190, 0, h_near, 247);
  gmsh::model::occ::addPoint(0.1411000614240227, 0.3501861850905401, 0, h_near, 248);
  gmsh::model::occ::addPoint(0.1546594447193058, 0.3484072716013699, 0, h_near, 249);
  gmsh::model::occ::addPoint(0.1683267575189567, 0.3464185151508629, 0, h_near, 250);
  gmsh::model::occ::addPoint(0.1820657524691879, 0.3442120317176733, 0, h_near, 251);
  gmsh::model::occ::addPoint(0.1958401822162117, 0.3417799372804556, 0, h_near, 252);
  gmsh::model::occ::addPoint(0.2096137994062407, 0.3391143478178640, 0, h_near, 253);
  gmsh::model::occ::addPoint(0.2233503566854871, 0.3362073793085530, 0, h_near, 254);
  gmsh::model::occ::addPoint(0.2370136067001635, 0.3330511477311768, 0, h_near, 255);
  gmsh::model::occ::addPoint(0.2505686505145814, 0.3296374883012909, 0, h_near, 256);
  gmsh::model::occ::addPoint(0.2640330089466618, 0.3259473215689799, 0, h_near, 257);
  gmsh::model::occ::addPoint(0.2774922979283390, 0.3219473895478331, 0, h_near, 258);
  gmsh::model::occ::addPoint(0.2910366843026320, 0.3176034866759807, 0, h_near, 259);
  gmsh::model::occ::addPoint(0.3047563349125602, 0.3128814073915527, 0, h_near, 260);
  gmsh::model::occ::addPoint(0.3187414166011429, 0.3077469461326795, 0, h_near, 261);
  gmsh::model::occ::addPoint(0.3330820962113991, 0.3021658973374913, 0, h_near, 262);
  gmsh::model::occ::addPoint(0.3478685405863481, 0.2961040554441182, 0, h_near, 263);
  gmsh::model::occ::addPoint(0.3631909165690094, 0.2895272148906906, 0, h_near, 264);
  gmsh::model::occ::addPoint(0.3791393910024019, 0.2824011701153385, 0, h_near, 265);
  gmsh::model::occ::addPoint(0.3958015194316177, 0.2746925523854085, 0, h_near, 266);
  gmsh::model::occ::addPoint(0.4132257846475764, 0.2663805144128162, 0, h_near, 267);
  gmsh::model::occ::addPoint(0.4314305911576642, 0.2574538479423018, 0, h_near, 268);
  gmsh::model::occ::addPoint(0.4504335697513631, 0.2479015926680023, 0, h_near, 269);
  gmsh::model::occ::addPoint(0.4702523512181548, 0.2377127882840553, 0, h_near, 270);
  gmsh::model::occ::addPoint(0.4909045663475211, 0.2268764744845982, 0, h_near, 271);
  gmsh::model::occ::addPoint(0.5124078459289438, 0.2153816909637684, 0, h_near, 272);
  gmsh::model::occ::addPoint(0.5347798207519048, 0.2032174774157031, 0, h_near, 273);
  gmsh::model::occ::addPoint(0.5580381216058857, 0.1903728735345399, 0, h_near, 274);
  gmsh::model::occ::addPoint(0.5822003792803687, 0.1768369190144160, 0, h_near, 275);
  gmsh::model::occ::addPoint(0.6072842245648352, 0.1625980868079891, 0, h_near, 276);
  gmsh::model::occ::addPoint(0.6333072882487674, 0.1476406612941676, 0, h_near, 277);
  gmsh::model::occ::addPoint(0.6602872011216466, 0.1319470495207083, 0, h_near, 278);
  gmsh::model::occ::addPoint(0.6882415939729553, 0.1154996496800321, 0, h_near, 279);
  gmsh::model::occ::addPoint(0.7171880975921747, 0.0982808599645601, 0, h_near, 280);
  gmsh::model::occ::addPoint(0.7471443427687871, 0.0802730785667134, 0, h_near, 281);
  gmsh::model::occ::addPoint(0.7781279602922739, 0.0614587036789129, 0, h_near, 282);
  gmsh::model::occ::addPoint(0.8101565809521170, 0.0418201334935799, 0, h_near, 283);
  gmsh::model::occ::addPoint(0.8432478355377984, 0.0213397662031352, 0, h_near, 284);
  gmsh::model::occ::addPoint(0.8774193548388000, 0.0000000000000000, 0, h_near, 285);

  // Surrounding boundary (points)

  gmsh::model::occ::addPoint(Lx, 0, 0, h_far, 290);
  gmsh::model::occ::addPoint(Lx, Ly, 0, h_far, 291);
  gmsh::model::occ::addPoint(L_duct, Ly, 0, h_far, 292);

  gmsh::model::occ::addLine(290, 285, 17);
  gmsh::model::occ::addLine(290, 291, 18);
  gmsh::model::occ::addLine(291, 292, 19);
  gmsh::model::occ::addLine(292, 1, 20);

  // Surrounding Cartesian PML
  if(Npml != 0) {
    gmsh::model::occ::addPoint(Lx + Npml * h_far, 0, 0, h_far, 293);
    gmsh::model::occ::addPoint(Lx + Npml * h_far, Ly + Npml * h_far, 0, h_far, 294);
    gmsh::model::occ::addPoint(L_duct - Npml * h_far, Ly + Npml * h_far, 0, h_far, 295);
    gmsh::model::occ::addPoint(L_duct - Npml * h_far, 1.5000000000000000, 0, h_near, 296);

    gmsh::model::occ::addLine(1, 296, 21);
    gmsh::model::occ::addLine(296, 295, 22);
    gmsh::model::occ::addLine(295, 294, 23);
    gmsh::model::occ::addLine(294, 293, 24);
    gmsh::model::occ::addLine(293, 290, 25);
  }
  // Nacelle surface (hard wall)
  std::vector< int > t;
  for(int i = 1; i <= 165; i++) t.push_back(i);
  gmsh::model::occ::addSpline(t, 1);

  // Nacelle surface (lined)
  t.clear();
  for(int i = 165; i <= 229; i++) t.push_back(i);
  gmsh::model::occ::addSpline(t, 2);

  // Nacelle surface (hard wall)
  gmsh::model::occ::addSpline({229, 230, 231, 232, 233}, 3);

  // Fan duct outer wall
  gmsh::model::occ::addLine(233, 234, 4);

  // Fan Face
  gmsh::model::occ::addLine(234, 235, 5);

  // Fan duct inner wall
  gmsh::model::occ::addLine(235, 236, 6);

  // Spinner surface
  t.clear();
  for(int i = 236; i <= 285; i++) t.push_back(i);
  gmsh::model::occ::addSpline(t, 7);

  if(Npml != 0) {
    //gmsh::model::occ::addCurveLoop({17,-7,-6,-5,-4,-3,-2,-1, -21, -22, -23, -24, -25}, 1);
    gmsh::model::occ::addCurveLoop({17, -7, -6, -5, -4, -3, -2, -1, -20, -19, -18}, 1);
    gmsh::model::occ::addCurveLoop({18, 19, 20, -21, -22, -23, -24, -25}, 2);
  }
  else {
    gmsh::model::occ::addCurveLoop({17, -7, -6, -5, -4, -3, -2, -1, -20, -19, -18}, 1);
  }

  if(Npml_active != 0) {
    gmsh::model::occ::addPoint(L_duct - Npml_active * h_far, 1.2000000000000000, 0, h_near, 500);
    gmsh::model::occ::addPoint(L_duct - Npml_active * h_far, 0.3586206896556000, 0, h_near, 501);
    gmsh::model::occ::addLine(235, 501, 502);
    gmsh::model::occ::addLine(501, 500, 503);
    gmsh::model::occ::addLine(500, 234, 504);
    gmsh::model::occ::addCurveLoop({503, 504, 5, 502}, 3);
  }

  gmsh::model::occ::addPlaneSurface({1}, 1);
  if(Npml != 0) {
    gmsh::model::occ::addPlaneSurface({2}, 2);
  }
  if(Npml_active != 0) {
    gmsh::model::occ::addPlaneSurface({3}, 3);
  }
  gmsh::model::occ::synchronize();

  // Physical groups
  // Volumes
  if(Npml != 0) {
    gmsh::model::addPhysicalGroup(2, {1, 2}, 1000);
    gmsh::model::setPhysicalName(2, 1000, "omega");
    gmsh::model::addPhysicalGroup(2, {1}, 2000);
    gmsh::model::setPhysicalName(2, 2000, "omega_phy");
    gmsh::model::addPhysicalGroup(2, {2}, 3000);
    gmsh::model::setPhysicalName(2, 3000, "omega_pml");
  }
  else {
    gmsh::model::addPhysicalGroup(2, {1}, 2000);
    gmsh::model::setPhysicalName(2, 2000, "omega");
  }

  if(Npml_active != 0) {
    gmsh::model::addPhysicalGroup(2, {3}, 4000);
    gmsh::model::setPhysicalName(2, 4000, "omega_pml_active");
    gmsh::model::addPhysicalGroup(1, {503}, 25000);
    gmsh::model::setPhysicalName(1, 25000, "FanFace_pml");
  }
  // Boundaries
  gmsh::model::addPhysicalGroup(1, {6, 7}, 6000);
  gmsh::model::setPhysicalName(1, 6000, "Spinner");


  gmsh::model::addPhysicalGroup(1, {2}, 2000);
  gmsh::model::setPhysicalName(1, 2000, "Liner");
  gmsh::model::addPhysicalGroup(1, {1, 3, 4}, 4000);
  gmsh::model::setPhysicalName(1, 4000, "NacelleSurface");

  gmsh::model::addPhysicalGroup(1, {5}, 5000);
  gmsh::model::setPhysicalName(1, 5000, "FanFace");

  gmsh::model::addPhysicalGroup(1, {17}, 17000);
  gmsh::model::setPhysicalName(1, 17000, "SymAxis");

  gmsh::model::addPhysicalGroup(1, {18}, 18000);
  gmsh::model::setPhysicalName(1, 18000, "ExtBoundary_Right");

  gmsh::model::addPhysicalGroup(1, {19}, 19000);
  gmsh::model::setPhysicalName(1, 19000, "ExtBoundary_Top");

  gmsh::model::addPhysicalGroup(1, {20}, 20000);
  gmsh::model::setPhysicalName(1, 20000, "ExtBoundary_Left");

  if(Npml != 0) {
    gmsh::model::addPhysicalGroup(1, {21, 22, 23, 24, 25}, 21000);
    gmsh::model::setPhysicalName(1, 21000, "ExtBoundary_Pml");
    gmsh::model::addPhysicalGroup(1, {22}, 22000);
    gmsh::model::setPhysicalName(1, 22000, "ExtBoundary_PmlLeft");
    gmsh::model::addPhysicalGroup(1, {24}, 24000);
    gmsh::model::setPhysicalName(1, 24000, "ExtBoundary_PmlRight");
  }

  gmsh::model::occ::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::model::mesh::setOrder(ElementOrder);
  gmsh::model::mesh::createEdges();
  gmsh::model::mesh::createFaces();

  gmsh::write("nacelle.msh");
}

void getInputMode(int m, int n, double w, double M, std::complex< double > &kzmn, double &krmn, ScalarFunction< std::complex< double > > &psi, ScalarFunction< std::complex< double > > &dr_psi, double Rin)
{

  ScalarFunction< std::complex< double > > Umn, dUmndr;
  ScalarFunction< std::complex< double > > Bmn, Bmn_var;
  std::complex< double > im(0., 1.);
  // get Laplace-Beltrami eigenvalue associated to the mode (m,n)
  double data[50][30];
  std::ifstream table("../Bessel_annular_zeros_m50_n30.txt"); // tabulated values for a given Rin, Rout !!
  if(Rin != 0.3586206896556000) { // tabulated values for a given Rin, Rout !!
    msg::error << "The Laplace-Beltrami eigenvalues must be recomputed !" << msg::endl;
  }
  if(!table) {
    std::cout << "Error, file couldn't be opened" << std::endl;
  }
  else {
    for(int row = 0; row < 50; row++) { // stop loops if nothing to read
      for(int column = 0; column < 30; column++) {
        table >> data[row][column];
        if(!table) {
          std::cout << "Error reading file for element " << row << "," << column << std::endl;
        }
      }
    }
  }
  krmn = data[m][n];

  double beta = sqrt(1 - M * M);
  double Argsqrt = w * w - beta * beta * krmn * krmn;
  if(Argsqrt >= 0) {
    kzmn = (1 / (beta * beta)) * (-M * w + sqrt(Argsqrt));
    msg::info << " - Input propagative mode " << msg::endl;
  }
  else {
    kzmn = (1 / (beta * beta)) * (-M * w - im * sqrt(abs(Argsqrt)));
    msg::info << " - Input evanescent mode " << msg::endl;
  }

  if(w > (beta * krmn) && w < krmn) {
    msg::info << " - Input inverse upstream mode ! " << msg::endl;
  }

  // compute Bmn
  ScalarFunction< std::complex< double > > dJmdr, dYmdr, dJmdr_var, dYmdr_var;
  ScalarFunction< std::complex< double > > kr = krmn * Rin;
  ScalarFunction< std::complex< double > > kr_var = krmn * y< std::complex< double > >();
  if(m == 0) {
    dJmdr = -cylBesselJ(1, kr);
    dYmdr = -cylNeumann(1, kr);
    dJmdr_var = -cylBesselJ(1, kr_var);
    dYmdr_var = -cylNeumann(1, kr_var);
  }
  else {
    dJmdr = 0.5 * (cylBesselJ(m - 1, kr) - cylBesselJ(m + 1, kr));
    dYmdr = 0.5 * (cylNeumann(m - 1, kr) - cylNeumann(m + 1, kr));
    dJmdr_var = 0.5 * (cylBesselJ(m - 1, kr_var) - cylBesselJ(m + 1, kr_var));
    dYmdr_var = 0.5 * (cylNeumann(m - 1, kr_var) - cylNeumann(m + 1, kr_var));
  }
  Bmn = -dJmdr / dYmdr;

  // compute Umn
  if(krmn == 0) {
    Umn = cylBesselJ(m, kr_var);
    dUmndr = krmn * dJmdr_var;
  }
  else {
    Umn = cylBesselJ(m, kr_var) + Bmn * cylNeumann(m, kr_var);
    dUmndr = krmn * (dJmdr_var + Bmn * dYmdr_var);
  }
  psi = exp(-im * kzmn * x< std::complex< double > >()) * Umn;
  dr_psi = exp(-im * kzmn * x< std::complex< double > >()) * dUmndr;
}
