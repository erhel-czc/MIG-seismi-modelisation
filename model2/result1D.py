import os
import pickle
import matplotlib.pyplot as plt
import numpy as np

class Result1D:
    "Results storage"

    def __init__(self, T, V, Nu, Phi, pd, pnd, pc, filename = ''):
        self.T=T
        self.V=V
        self.Nu=Nu
        self.Phi=Phi
        self.pd=pd
        self.pnd=pnd
        self.pc=pc

        if filename=='' or filename is None:
            # Build a robust default filename using available scalar attributes.
            try:
                # pnd.a may be an array: take its mean as representative
                a_val = float(np.mean(pnd.a)) if hasattr(pnd, 'a') else 0.0
            except Exception:
                a_val = 0.0

            eta_val = float(getattr(pnd, 'eta', 0.0))
            h_val = float(getattr(pnd, 'h', 0.0))

            self.filename = f"a{a_val:.2e}_eta{eta_val:.2e}_h{h_val:.2e}.pkl"
        else:
            self.filename = filename

    def save_results(self, folder_name=''):
        # Save results under the package-relative `model2/Results1D/<folder_name>/`
        base_dir = os.path.dirname(__file__)
        results_dir = os.path.join(base_dir, 'Results1D')

        if folder_name:
            results_dir = os.path.join(results_dir, folder_name)

        os.makedirs(results_dir, exist_ok=True)

        outpath = os.path.join(results_dir, self.filename)
        
        with open(outpath, 'wb') as f:
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

        """
        # The following code was found on internet sources, in order to bypass
        # issues with unpickling classes.

        # When objects were pickled from a script run as __main__,
        # pickle stores the classes with module name '__main__'.
        # When loading from another script (like this one) we need
        # to make sure those class names are available on the
        # current '__main__' module so unpickling can find them.
        import sys
        
        import 1D_CFAULT_RNS as prmod

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
                    pass"""

        with open(path, 'rb') as f:
            data = pickle.load(f)

        return data

    def theta_evolution(self, save = False, path = '', name=''):
        """Plot state variable evolution"""
        plt.figure('Theta evolution')

        plt.plot(self.T,self.Nu,'-k')
        plt.xlabel('Time (ND)')
        plt.ylabel('Log of state variable (ND)')
        plt.title(r'$\theta$ evolution')
        plt.grid()

        if save:
            if name != '':
                plt.savefig(f"{path}/theta_evolution_{name}.png", dpi=300)
            else:
                plt.savefig(f"{path}/theta_evolution.png", dpi=300)