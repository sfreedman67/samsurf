from bowman import triangulation as trin


class TestAlgo:
    def test_generators_veech(self):
        x = trin.Triangulation.regular_octagon()
        f = x.generators_veech
        assert False
