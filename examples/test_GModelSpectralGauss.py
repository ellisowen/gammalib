import gammalib

def test_eval():
    model = gammalib.GModelSpectralGauss()
    energy = gammalib.GEnergy(42, "MeV")
    time = gammalib.GTime()
    print(model.eval(energy, time))
    
def test_mc():
    model = gammalib.GModelSpectralGauss()
    emin = gammalib.GEnergy(1, "MeV")
    emax = gammalib.GEnergy(3, "MeV")
    time = gammalib.GTime()
    ran = gammalib.GRan()
    print(model.mc(emin, emax, time, ran))
    
def test_write():
    gammalib.GModelSpectralGauss.write
    pass

def test_flux():
    gammalib.GModelSpectralGauss.flux
    pass


if __name__ == '__main__':
    test_eval()
    test_mc()
    