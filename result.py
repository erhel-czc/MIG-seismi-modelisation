import pickle
import matplotlib.pyplot as plt
import numpy as np

class Result:
    "Results storage"

    def __init__(self, T, V, Vpoint, Nu, Phi, Phipoint, Tau, Sigma_n, pd, pnd, pc, filename = ''):
        self.T=T
        self.V=V
        self.Vpoint=Vpoint
        self.Nu=Nu
        self.Phi=Phi
        self.Phipoint=Phipoint
        self.Tau=Tau
        self.Sigma_n=Sigma_n
        self.pd=pd
        self.pnd=pnd
        self.pc=pc

        if filename=='':
            if Tau is not [] and Sigma_n is not []:
                self.filename= f"{pnd.a}_{pnd.eta}_{pnd.k}_{Sigma_n[0]}.pkl"

            else:
                self.filename= f"{pnd.a}_{pnd.eta}_{pnd.k}.pkl"
        else:
            self.filename=filename

    def save_results(self, folder_name=''):
        with open(f"Results/{folder_name}/{self.filename}", 'wb') as f:
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

    def slip_rate_evolution(self, save=False, path=''):
        """
        Plot the slip rate evolution of the speed.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        """
        plt.figure('Slip rate evolution')

        plt.plot(self.T, self.Phi, '-k')
        plt.xlabel('Time (ND)')
        plt.ylabel('Log Slip rate (ND)')
        plt.title(r'Slip rate evolution ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
        plt.grid()

        # plt.xlim([0, 100])
        #plt.savefig('images/modif_parameters/slip_rate_k%.2f_a%.2f.pdf' % (self.pnd.k, self.pnd.a))
        if save:
            plt.savefig(f'{path}/slip_rate_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')

    def phase_portrait(self, save=False, path=''):
        """
        Plot the phase portrait of the speed.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        """
        plt.figure('Phase portrait')

        sc = plt.scatter(self.V[1:], self.Vpoint, c=self.T[1:], cmap='viridis', marker='+')
        plt.plot(self.V[1:], self.Vpoint, '-k', alpha=0.3)
        plt.xlabel('Speed (ND)')
        plt.ylabel('Acceleration (ND)')
        plt.title(r'Phase portrait ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
        plt.colorbar(sc, label='Time progression')
        plt.grid()

        if save:
            plt.savefig(f'{path}/phase_portrait_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')

    def tau_evolution(self, save=False, path=''):
        """
        Plot the tangent stress evolution.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        """
        plt.figure('Shear stress evolution')

        plt.plot(self.T, self.Tau, '-k')
        plt.xlabel('Time (ND)')
        plt.ylabel('Shear stress (ND)')
        plt.title(r'Shear stress evolution ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
        plt.grid()

        if save:
            plt.savefig(f'{path}/shear_stress_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')

    def sigma_evolution(self, save=False, path=''):
        """
        Plot the normal stress evolution.

        Parameters
        ----------
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        """
        plt.figure('Normal stress evolution')

        plt.plot(self.T, self.Sigma_n, '-k')
        plt.xlabel('Time (ND)')
        plt.ylabel('Normal stress (ND)')
        plt.title(r'Normal stress evolution ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
        plt.grid()

        if save:
            plt.savefig(f'{path}/normal_stress_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')
