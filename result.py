import pickle
import matplotlib.pyplot as plt
import numpy as np

class Result:
    "Results storage"

    def __init__(self, T, V, Vpoint, Nu, Phi, Phipoint, Tau, Sigma_n, pd, pnd, pc, P=None, filename = ''):
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
        self.P = P

        if filename=='':
            if Sigma_n is not None:
                self.filename = f"{pnd.a:.2e}_{pnd.eta:.2e}_{pnd.k:.2e}_{Sigma_n[0]:.2e}.pkl"
            else:
                self.filename= f"{pnd.a:.2e}_{pnd.eta:.2e}_{pnd.k:.2e}.pkl"
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
    
    def magnitude(self, print_result=True):
        try :
            threshold_velocity = 1/self.pd.V_p  # Define a threshold velocity to identify slip events

            T=self.T
            V=self.V
            deplacement=0

            if max(V) > threshold_velocity:
                while len(T)>0 and V[0] < threshold_velocity:
                    T = T[1:]
                    V = V[1:]

                while len(T)>0 and V[0] > threshold_velocity:
                    deplacement += V[0] * (T[1] - T[0])
                    T = T[1:]
                    V = V[1:]
                
                deplacement *= self.pd.dc
                
                if print_result:
                    print(f"Slip event detected. Total slip displacement: {deplacement}.")   

                def compute_magnitude(deplacement, pd):
                    """
                    Calculate the moment magnitude based on slip displacement and physical parameters.

                    Parameters:
                    deplacement (float): Slip displacement in meters.
                    pd (PhysicalData): An instance of the PhysicalData class containing physical parameters.

                    Returns:
                    float: Calculated moment magnitude.
                    """
                    mu = pd.shear
                    L = pd.lenght_fault

                    # Calculate seismic moment (M0)
                    M0 = mu * L**2 * deplacement  # in Newton-meters

                    # Calculate moment magnitude (Mw)
                    Mw = (2/3) * np.log10(M0) - 6.0

                    return Mw

                if print_result :print(f"Calculated moment magnitude: {compute_magnitude(deplacement, self.pd)} for {self.filename}.")

                return compute_magnitude(deplacement, self.pd)
            else:
                if print_result : print(f"No slip event detected for {self.filename}.")
    
        except :
            print(f"Magnitude calculation failed for {self.filename}.")
    
    def cycles(self, print_result=True):
        T=self.T
        Phi=self.Phi
        maxima=[]
        periods=[]

        while len(T)>1:
            while len(T)>1 and Phi[0] < Phi[1]:
                T = T[1:]
                Phi = Phi[1:]
            maxima.append(T[0])
            T = T[1:]
            Phi = Phi[1:]
            while len(T)>1 and Phi[0] > Phi[1]:
                T = T[1:]
                Phi = Phi[1:]

        for i in range(1, len(maxima)):
            periods = np.append(periods, [maxima[i]-maxima[i-1]])
        
        if print_result :
            if len(periods)>0:
                print(f"Cycles: {periods*self.pd.dc/self.pd.V_p} in s for {self.filename}.")
            else:
                print(f"No cycle detected for {self.filename}.")


    def slip_rate_evolution(self, save=False, path='', name=''):
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
            if name == '':
                plt.savefig(f'{path}/slip_rate_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')
            else:
                plt.savefig(f'{path}/slip_rate_{name}.pdf')

    def phase_portrait(self, save=False, path='', name=''):
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
            if name == '':
                plt.savefig(f'{path}/phase_portrait_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')
            else:
                plt.savefig(f'{path}/phase_portrait_{name}.pdf')

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
        plt.figure('Tau stress evolution')

        plt.plot(self.T*self.pd.dc/self.pd.V_p, self.Tau*self.pd.sigma_n0, '-k')
        plt.xlabel('Time (s)')
        plt.ylabel('Shear stress (Pa)')
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

        plt.plot(self.T*self.pd.dc/self.pd.V_p, self.pd.sigma_n0*self.Sigma_n, '-k')
        plt.xlabel('Time (s)')
        plt.ylabel('Normal stress (Pa)')
        plt.title(r'Normal stress evolution ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
        plt.grid()

        if save:
            plt.savefig(f'{path}/normal_stress_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')
        
    def theta_evolution(self, save=False, path='', name=''):
            """
            Plot the theta evolution.

            Parameters
            ----------
            save : bool
                If True, save the figure in the given path.
            path : str
                The path where the figure will be saved.
            """
            plt.figure('Log Theta evolution')

            plt.plot(self.T, self.Nu, '-k')
            plt.xlabel('Time (ND)')
            plt.ylabel('Log Theta (ND)')
            plt.title(r'Log Theta evolution ($\kappa$=%.2f, $\alpha$=%.2f)' % (self.pnd.k, self.pnd.a))
            plt.grid()

            if save:
                if name == '':
                    plt.savefig(f'{path}/normal_stress_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')
                else:
                    plt.savefig(f'{path}/theta_evol_{name}.pdf')

    def pressure_evolution(self, save=False, path='', name=''):

        plt.figure('Pressure evolution')

        plt.plot(self.T*self.pd.dc/self.pd.V_p, self.P, '-k') # type: ignore
        plt.xlabel('Time (s)')
        plt.ylabel('Pressure (Pa)')
        plt.title(r'Pressure evolution')
        plt.grid()

        if save:
            if name == '':
                plt.savefig(f'{path}/pressure_evol_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')
            else:
                plt.savefig(f'{path}/pressure_evol_{name}.pdf')

    def fft_cycle_period(self, output=False, threshold_amplitude=0.05, save=False, path='', name=''):
        """
        Plot the FFT of the slip rate signal to determine the cycle period.

        Parameters
        ----------
        output : bool
            If True, print the dominant cycle period, and show the figure.
        save : bool
            If True, save the figure in the given path.
        path : str
            The path where the figure will be saved.
        name : str
            The name to use when saving the figure.
        """
        signal = self.V
        time = self.T

        sampling_rate = 1 / np.mean(np.diff(time))
        fft = np.fft.rfft(signal)
        freqs = np.fft.rfftfreq(len(signal), d=1/sampling_rate)
        
        period = 1 / freqs[np.argmax(np.abs(fft[1:])) + 1]

        max_amp = np.max(np.abs(fft))
        threshold = threshold_amplitude * max_amp
        mask = np.abs(fft) >= threshold
        fft = fft[mask]
        freqs = freqs[mask]

        
        if output:
            print(f"Dominant cycle period: {period} (ND) for k={self.pnd.k}, a={self.pnd.a}.")
            plt.figure('FFT of v(t)')
            plt.plot(freqs, np.abs(fft), '-k')
            plt.xlabel('Frequency (Hz)')
            plt.ylabel('Amplitude')
            plt.title(r'FFT of v(t), cycle period T=%.2f, ($\kappa$=%.2f, $\alpha$=%.2f)' % (period, self.pnd.k, self.pnd.a))
            plt.grid()
            plt.show()

        if save:
            if name == '':
                plt.savefig(f'{path}/fft_cycle_period_k{round(self.pnd.k,2)}_a{round(self.pnd.a,2)}.pdf')
            else:
                plt.savefig(f'{path}/fft_cycle_period_{name}.pdf')

        return period
