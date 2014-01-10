import gammalib

def test_a():
    model = gammalib.GModelSpectralGauss()
    energy = gammalib.GEnergy(1, "TeV")
    time = gammalib.GTime()
    print(model.eval(energy, time))


if __name__ == '__main__':
    test_a()