// This custom random function is a macro to replace a random call from 0 to X but with floats
#define UNIFORM_RAND_0_TO_X(X) ((double)rand() / RAND_MAX * X)

// NOTE: before, it was calling rand() * 2 - 1, which actually gave enormous numbers, not floats
// I assumed that these were supposed to be random float calls between 0 and 2 - 1 so that we get numbers
// uniformly between -1 and 1
double amplitude = 0.000001;                               //10^(-6)
var[VAR::U_ALPHA] = 1 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1); // lapse
var[VAR::U_CHI] = 1 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // chi
var[VAR::U_K] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);          // trace K
var[VAR::U_CAP_GT0] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);    // Gt0
var[VAR::U_CAP_GT1] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);    // Gt1
// Gt2
var[VAR::U_CAP_GT2] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);
var[VAR::U_BETA0] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);       //2. * a * m * o14 * o148 * o21 * o28 * yy -> 0
var[VAR::U_BETA1] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);       // -2. * a * m * o14 * o148 * o21 * o28 * xx -> 0
var[VAR::U_BETA2] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);       // shift2
var[VAR::U_B0] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);          // gaugeB0
var[VAR::U_B1] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);          // gaugeB1
var[VAR::U_B2] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);          // gaugeB2
var[VAR::U_GT00] = 1 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        //gt11
var[VAR::U_GT01] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // gt12
var[VAR::U_GT02] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // gt13
var[VAR::U_GT11] = 1 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // gt22
var[VAR::U_GT12] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // gt23
var[VAR::U_GT22] = 1 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // gt33
var[VAR::U_AT00] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // At11
var[VAR::U_AT01] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // At12
var[VAR::U_AT02] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // At13
var[VAR::U_AT11] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // At22
var[VAR::U_AT12] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // At23
var[VAR::U_AT22] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);        // At33
var[VAR::U_DILATONPHI] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);  // dilatonPhi
var[VAR::U_KAPPA] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);       // -2. * a * m * o174 * o71 * o9 * zz -> 0
var[VAR::U_CAPITALPI] = 0. + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);  // capitalPi
var[VAR::U_CAPITALXI] = 0. + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);  // capitalXi
var[VAR::U_DAMPINGPSI] = 0. + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1); // dampingPsi
var[VAR::U_DAMPINGPHI] = 0. + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1); // dampingPhi
var[VAR::U_PERPE0] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);      // perpE0
var[VAR::U_PERPE1] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);      // perpE1
var[VAR::U_PERPE2] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);      // perpE2
var[VAR::U_PERPB0] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);      //-0.7071067811865475 * a * m * o148 * o179 * o180 * o188 * o196 * xx * zz -> 0
var[VAR::U_PERPB1] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);      //-0.7071067811865475 * a * m * o148 * o179 * o180 * o188 * o196 * yy * zz -> 0
var[VAR::U_PERPB2] = 0 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);      //-0.7071067811865475 * a * m * o148 * o179 * o180 * o188 * (2. * o14 * o185 * o6 - o18 * o183 * o61) -> 0
