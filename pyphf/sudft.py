import util


class SUDFT():
    def __init__(self, guesshf):
        suhf = util.SUHF(guesshf)
        self.suhf = suhf

    def kernel(self):
        self.suhf.kernel()
