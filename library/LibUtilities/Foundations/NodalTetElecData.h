
namespace Nektar
{
namespace LibUtilities
{
static const unsigned int perm4_3d[4][4] = {
    {0, 1, 2, 3},
    {3, 0, 1, 2},
    {2, 3, 0, 1},
    {1, 2, 3, 0}}; // Works for abbb or aaab
static const unsigned int perm6_3d[6][4] = {
    {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 2, 3, 1},
    {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 3, 0, 1}}; // Works for aabb

static const unsigned int perm12A_3d[12][4] = {
    {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1},
    {0, 3, 1, 2}, {0, 3, 2, 1}, {2, 0, 1, 3}, {2, 0, 3, 1},
    {2, 3, 0, 1}, {3, 0, 1, 2}, {3, 0, 2, 1}, {3, 2, 0, 1}}; // Works for aabc
static const unsigned int perm12B_3d[12][4] = {
    {0, 1, 2, 3}, {0, 2, 1, 3}, {0, 2, 3, 1}, {1, 0, 2, 3},
    {1, 2, 0, 3}, {1, 2, 3, 0}, {2, 0, 1, 3}, {2, 0, 3, 1},
    {2, 1, 0, 3}, {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0}}; // Works for abcc
static const unsigned int perm12C_3d[12][4] = {
    {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 3, 1, 2}, {1, 0, 2, 3},
    {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0}, {1, 3, 0, 2},
    {1, 3, 2, 0}, {3, 0, 1, 2}, {3, 1, 0, 2}, {3, 1, 2, 0}}; // Works for abbc

static const unsigned int perm24_3d[24][4] = {
    {0, 1, 2, 3}, {0, 1, 3, 2}, {0, 2, 1, 3}, {0, 2, 3, 1}, {0, 3, 1, 2},
    {0, 3, 2, 1}, {1, 0, 2, 3}, {1, 0, 3, 2}, {1, 2, 0, 3}, {1, 2, 3, 0},
    {1, 3, 0, 2}, {1, 3, 2, 0}, {2, 0, 1, 3}, {2, 0, 3, 1}, {2, 1, 0, 3},
    {2, 1, 3, 0}, {2, 3, 0, 1}, {2, 3, 1, 0}, {3, 0, 1, 2}, {3, 0, 2, 1},
    {3, 1, 0, 2}, {3, 1, 2, 0}, {3, 2, 0, 1}, {3, 2, 1, 0}}; // Works for abcd

const unsigned int NodalTetElecAvailable                          = 10;
static const unsigned int NodalTetElecNPTS[NodalTetElecAvailable] = {
    1, 2, 3, 5, 6, 9, 11, 15, 18, 23};
static const NekDouble NodalTetElecData[][9] = {
    // %%% n_1    n_4    n_6   n_12    n_24          l_1             l_2 l_3 l_4
    // 1 1 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    // 2 2 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {0, 0, 1, 0, 0, 0.5000000000, 0.5000000000, 0.0000000000, 0.0000000000},
    // 3 3 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {0, 1, 0, 0, 0, 0.3333333333, 0.3333333333, 0.3333333333, 0.0000000000},
    {0, 0, 0, 2, 0, 0.7236067977, 0.2763932023, 0.0000000000, 0.0000000000},
    // 4 5 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {1, 0, 0, 0, 0, 0.2500000000, 0.2500000000, 0.2500000000, 0.2500000000},
    {0, 0, 1, 0, 0, 0.5000000000, 0.5000000000, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.8273268354, 0.1726731646, 0.0000000000, 0.0000000000},
    {0, 0, 0, 1, 0, 0.2371200168, 0.2371200168, 0.5257599664, 0.0000000000},
    // 5 6 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {0, 1, 0, 0, 0, 0.1834903473, 0.1834903473, 0.1834903473, 0.4495289581},
    {0, 0, 0, 2, 0, 0.8825276620, 0.1174723380, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.6426157582, 0.3573842418, 0.0000000000, 0.0000000000},
    {0, 0, 0, 1, 0, 0.1575181512, 0.1575181512, 0.6849636976, 0.0000000000},
    {0, 0, 0, 1, 0, 0.4105151510, 0.4105151510, 0.1789696980, 0.0000000000},
    // 6 9 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {0, 1, 0, 0, 0, 0.3333333333, 0.3333333333, 0.3333333333, 0.0000000000},
    {0, 1, 0, 0, 0, 0.1402705801, 0.1402705801, 0.1402705801, 0.5791882597},
    {0, 0, 1, 0, 0, 0.5000000000, 0.5000000000, 0.0000000000, 0.0000000000},
    {0, 0, 1, 0, 0, 0.3542052583, 0.3542052583, 0.1457947417, 0.1457947417},
    {0, 0, 0, 2, 0, 0.9151119481, 0.0848880519, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.7344243967, 0.2655756033, 0.0000000000, 0.0000000000},
    {0, 0, 0, 1, 0, 0.1061169285, 0.1061169285, 0.7877661430, 0.0000000000},
    {0, 0, 0, 0, 1, 0.3097982151, 0.5569099204, 0.1332918645, 0.0000000000},
    // 7 11 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {0, 1, 0, 0, 0, 0.1144606542, 0.1144606542, 0.1144606542, 0.6566180374},
    {0, 1, 0, 0, 0, 0.2917002822, 0.2917002822, 0.2917002822, 0.1248991534},
    {0, 0, 0, 2, 0, 0.9358700743, 0.0641299257, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.7958500907, 0.2041499093, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.6046496090, 0.3953503910, 0.0000000000, 0.0000000000},
    {0, 0, 0, 1, 0, 0.0660520784, 0.0660520784, 0.8678958432, 0.0000000000},
    {0, 0, 0, 1, 0, 0.4477725053, 0.4477725053, 0.1044549894, 0.0000000000},
    {0, 0, 0, 1, 0, 0.2604038024, 0.2604038024, 0.4791923952, 0.0000000000},
    {0, 0, 0, 1, 0, 0.1208429970, 0.1208429970, 0.4770203357, 0.2812936703},
    {0, 0, 0, 0, 1, 0.2325524777, 0.6759625951, 0.0914849272, 0.0000000000},
    // 8 15 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {1, 0, 0, 0, 0, 0.2500000000, 0.2500000000, 0.2500000000, 0.2500000000},
    {0, 1, 0, 0, 0, 0.0991203900, 0.0991203900, 0.0991203900, 0.7026388300},
    {0, 0, 1, 0, 0, 0.5000000000, 0.5000000000, 0.0000000000, 0.0000000000},
    {0, 0, 1, 0, 0, 0.3920531037, 0.3920531037, 0.1079468963, 0.1079468963},
    {0, 0, 0, 2, 0, 0.9498789977, 0.0501210023, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.8385931398, 0.1614068602, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.6815587319, 0.3184412681, 0.0000000000, 0.0000000000},
    {0, 0, 0, 1, 0, 0.0660520784, 0.0660520784, 0.8678958432, 0.0000000000},
    {0, 0, 0, 1, 0, 0.2033467796, 0.2033467796, 0.5933064408, 0.0000000000},
    {0, 0, 0, 1, 0, 0.3905496216, 0.3905496216, 0.2189007568, 0.0000000000},
    {0, 0, 0, 1, 0, 0.1047451941, 0.1047451941, 0.5581946462, 0.2323149656},
    {0, 0, 0, 1, 0, 0.2419418605, 0.2419418605, 0.4062097450, 0.1099065340},
    {0, 0, 0, 0, 1, 0.3617970895, 0.5541643672, 0.0840385433, 0.0000000000},
    {0, 0, 0, 0, 1, 0.1801396087, 0.7519065566, 0.0679538347, 0.0000000000},
    // 9 18 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {0, 1, 0, 0, 0, 0.3333333333, 0.3333333333, 0.3333333333, 0.0000000000},
    {0, 1, 0, 0, 0, 0.0823287303, 0.0823287303, 0.0823287303, 0.7530138091},
    {0, 1, 0, 0, 0, 0.2123055477, 0.2123055477, 0.2123055477, 0.3630833569},
    {0, 0, 0, 2, 0, 0.9597669541, 0.0402330459, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.8693869326, 0.1306130674, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.7389624749, 0.2610375251, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.5826394788, 0.4173605212, 0.0000000000, 0.0000000000},
    {0, 0, 0, 1, 0, 0.0355775717, 0.0355775717, 0.9288448566, 0.0000000000},
    {0, 0, 0, 1, 0, 0.4640303025, 0.4640303025, 0.0719393950, 0.0000000000},
    {0, 0, 0, 1, 0, 0.1633923069, 0.1633923069, 0.6732153862, 0.0000000000},
    {0, 0, 0, 1, 0, 0.0873980781, 0.0873980781, 0.6297057875, 0.1954980564},
    {0, 0, 0, 1, 0, 0.0916714679, 0.0916714679, 0.4819523024, 0.3347047619},
    {0, 0, 0, 1, 0, 0.2040338880, 0.2040338880, 0.4996292993, 0.0923029247},
    {0, 0, 0, 1, 0, 0.3483881173, 0.3483881173, 0.2075502723, 0.0956734931},
    {0, 0, 0, 0, 1, 0.2966333890, 0.6349633653, 0.0684032457, 0.0000000000},
    {0, 0, 0, 0, 0, 0.1439089974, 0.8031490682, 0.0529419344, 0.0000000000},
    {0, 0, 0, 0, 0, 0.3225890045, 0.4968009397, 0.1806100558, 0.0000000000},
    // 10 23 %%% Order / Number of Points
    {0, 1, 0, 0, 0, 1.0000000000, 0.0000000000, 0.0000000000, 0.0000000000},
    {0, 1, 0, 0, 0, 0.0678316144, 0.0678316144, 0.0678316144, 0.7965051568},
    {0, 1, 0, 0, 0, 0.1805746957, 0.1805746957, 0.1805746957, 0.4582759129},
    {0, 1, 0, 0, 0, 0.3051527124, 0.3051527124, 0.3051527124, 0.0845418628},
    {0, 0, 1, 0, 0, 0.5000000000, 0.5000000000, 0.0000000000, 0.0000000000},
    {0, 0, 1, 0, 0, 0.3164336236, 0.3164336236, 0.1835663764, 0.1835663764},
    {0, 0, 1, 0, 0, 0.4219543801, 0.4219543801, 0.0780456199, 0.0780456199},
    {0, 0, 0, 2, 0, 0.9670007152, 0.0329992848, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.8922417368, 0.1077582632, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.7826176635, 0.2173823365, 0.0000000000, 0.0000000000},
    {0, 0, 0, 2, 0, 0.6478790678, 0.3521209322, 0.0000000000, 0.0000000000},
    {0, 0, 0, 1, 0, 0.0265250690, 0.0265250690, 0.9469498620, 0.0000000000},
    {0, 0, 0, 1, 0, 0.1330857076, 0.1330857076, 0.7338285848, 0.0000000000},
    {0, 0, 0, 1, 0, 0.4232062312, 0.4232062312, 0.1535875376, 0.0000000000},
    {0, 0, 0, 1, 0, 0.2833924371, 0.2833924371, 0.4332151258, 0.0000000000},
    {0, 0, 0, 1, 0, 0.1734555313, 0.1734555313, 0.5762731177, 0.0768158196},
    {0, 0, 0, 1, 0, 0.0724033935, 0.0724033935, 0.6893564961, 0.1658367169},
    {0, 0, 0, 1, 0, 0.0768451848, 0.0768451848, 0.5573732958, 0.2889363346},
    {0, 0, 0, 0, 1, 0.3934913008, 0.5472380443, 0.0592706549, 0.0000000000},
    {0, 0, 0, 0, 1, 0.2462883939, 0.6991456238, 0.0545659823, 0.0000000000},
    {0, 0, 0, 0, 1, 0.1163195334, 0.8427538829, 0.0409265838, 0.0000000000},
    {0, 0, 0, 0, 1, 0.2707097521, 0.5811217960, 0.1481684519, 0.0000000000},
    {0, 0, 0, 0, 1, 0.3019928872, 0.4393774966, 0.1776946096, 0.0809350066}};

} // namespace LibUtilities
} // namespace Nektar
