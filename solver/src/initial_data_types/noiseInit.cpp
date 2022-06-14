// This custom random function is a macro to replace a random call from 0 to X
// but with floats
#define UNIFORM_RAND_0_TO_X(X) ((double)rand() / RAND_MAX * X)

// NOTE: before, it was calling rand() * 2 - 1, which actually gave enormous
// numbers, not floats I assumed that these were supposed to be random float
// calls between 0 and 2 - 1 so that we get numbers uniformly between -1 and 1
double amplitude = 0.000001;                                         // 10^(-6)
var[VAR::U_ALPHA] = 1 + amplitude * (UNIFORM_RAND_0_TO_X(2) - 1);    // lapse
// ........ the rest
