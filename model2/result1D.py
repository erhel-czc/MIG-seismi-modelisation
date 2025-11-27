import pickle
import matplotlib.pyplot as plt
import numpy as np

class Result:
    "Results storage"

    def __init__(self, T, V, Nu, Phi, pd, pnd, pc, filename = ''):
        self.T=T
        self.V=V
        self.Nu=Nu
        self.Phi=Phi
        self.pd=pd
        self.pnd=pnd
        self.pc=pc
        if filename=='':
            self.filename= f"{pnd.a:.2e}_{pnd.eta:.2e}_{pnd.k:.2e}.pkl"
        else:
            self.filename=filename

    def save_results(self, folder_name=''):
        # print folders
        import os
        print(os.listdir())
        with open(f"model2/Results1D/{folder_name}/{self.filename}", 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def load_results(path=''):
        """
        Load previously saved results from a pickle file.

        Parameters
        ----------
        filename : str
            The name of the pickle file to load.

        Returns
        -------
        data : Result
            The loaded results.
        """

        # The following code was found on internet sources, in order to bypass
        # issues with unpickling classes.

        # When objects were pickled from a script run as __main__,
        # pickle stores the classes with module name '__main__'.
        # When loading from another script (like this one) we need
        # to make sure those class names are available on the
        # current '__main__' module so unpickling can find them.
        import sys
        
        import PR_QDYN_RNS as prmod

        main_mod = sys.modules.get('__main__')

        # Export any classes defined in PR_QDYN_RNS to __main__ so that
        # objects pickled when that script was run as __main__ can be
        # located during unpickling here.
        for name, obj in vars(prmod).items():
            if isinstance(obj, type):
                try:
                    setattr(main_mod, name, obj)
                except Exception:
                    # If for some reason we can't set an attribute,
                    # ignore and continue with other classes.
                    pass

        with open(path, 'rb') as f:
            data = pickle.load(f)

        return data