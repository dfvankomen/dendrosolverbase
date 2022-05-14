/**
 * This file is just included *as is* by the compiler to calculate the KO derivatives necessary for KO disappation
 *
 *
 */

/*[[[cog
import cog
import sys
import importlib.util
import dendrosym

# the following lines will import any module directly from
spec = importlib.util.spec_from_file_location("dendroconf", CONFIG_FILE_PATH)
dendroconf = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = dendroconf
spec.loader.exec_module(dendroconf)

cog.outl('//// KO GRADIENT CALCULATION')
cog.outl(dendroconf.dendroConfigs.generate_ko_derivs("evolution"))

]]]*/

//[[[end]]]