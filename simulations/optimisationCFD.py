import numpy as np
from CAD.regenCooling import regenNozzle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class optimisationCFD:
    def __init__(self, config, fuel, injector, nozzle, oxidizer): #dodati jos potrebnih klasa koje ce se ubaciti ovdje
        self.config = config
        self.fuel = fuel
        self.injector = injector
        self.nozzle = nozzle
        self.oxidizer = oxidizer
        self.dynamic_viscosity = fuel.get_dynamic_viscosity()
        self.density = fuel.get_density_fuel()
        self.epsilon = 5e-5          # m - roughness of the channel walls recalibrate for specific 3d printer
        self.chanel_length = 3 * config.combustion_chamber_height/1000  # m - length of the cooling channel
    
        self.n_max = 1
        self.diamteres = []
        while True:
            diam_max = self.find_minimum_diameter(self.n_max)
            if diam_max is None:
                #print(f"Physical/design limit reached at {self.n_max} channels.")
                self.n_max -= 1  # Decrease the number of channels if no suitable diameter is found
                break
            # print(f"maximum number of channels: {self.n_max} for channel diameter: {diam_max * 1000:.2f} mm")
            self.diamteres.append(diam_max)
            self.n_max += 1
        #print(self.n_max)
        resolution = self.n_max

        self.A1 = np.pi * (np.min(self.diamteres))**2 / 4  # m^2 za jedan kanal
        self.Amax = np.pi * (np.max(self.diamteres)*10)**2 / 4  # m^2 za jedan kanal

        #setup of 3d space for optimisation
        self.metal_thickness = np.linspace(1e-3, config.ALU_THICKNESS * 10, resolution)
        self.number_of_channels = np.linspace(1, self.n_max, resolution)
        self.povrsina_kanala = np.linspace(self.A1, self.Amax, resolution) 

        #dobijen 3d prosto, potrebno ga odrezati da se dobiju fizikalni rezultati (mozda)
        T, N, A = np.meshgrid(self.metal_thickness, self.number_of_channels, self.povrsina_kanala, indexing='ij')
        self.domain_points = np.column_stack((T.flatten(), N.flatten(), A.flatten()))
       
      
        # prije toga potrebno definirati minimalne rubne uvjete za optimizaciju kako bi se mogli dobiti testirati fizikalno dobri dizajnovi

        # generiranje dizajnova mlaznica
        self.regenNozzleDesign = regenNozzle(
            combustion_chamber_diameter=self.nozzle.calculate_combustion_chamber_diameter() + 2 * config.CHAMBERSAFE_THICKNESS/1000,
            vanjski_promjer_ulaz=config.OXIDIZER_TANK_DIAMETER,
            promjer_grla_mlaznice=self.nozzle.circle_diameter(self.nozzle.calculate_throat_area()) * 1000,
            tlak_komore=config.COMBUSTION_CHAMBER_PRESSURE,
            tlak_fluida=self.oxidizer.get_pressure_liquid(),
            cvrstoca_materijala=config.TANK_MATERIAL_STRENGTH,
            faktor_sigurnosti=config.SAFETY_FACTOR_CHAMBER,
            d_hidraulicke_cijevi=self.povrsina_kanala[10], # da se dobije Auk potrebno self.povrsina_kanala * self.number_of_channels
            maseni_protok_IPA=self.injector.fuel_massflowrate,
            gustoca_IPA=785,
            alu_thickness=self.metal_thickness[3],
            n=self.number_of_channels[11],
            k=1
        )

        # # Replace 'myscript.sh' with your script path
        # result = subprocess.run(['sbatch', 'myscript.sh'], capture_output=True, text=True)

        # print("STDOUT:", result.stdout)
        # print("STDERR:", result.stderr)

        #for i in range(9):
        
            #algo below must find 9 points in 3d space so design sphere can be visualised

            # 1. func to use randomised point and extra 3 points in its proximity in 3D space
            # 1.1. check if dimensions will fit for production, if not skip this area and select new random point (maybe this point is not needed)

            #for i in range(until it reaches desired temperature from cfd simulation):

                # 2. func to start cfd simulationto get temperatures from this 4 points in 3D space

                # 3. fucntion to calculate vector of temperatures from 4 points in 3D space

                # 4. fucntion to select new points in 3D space based on temperature vector

        #after 9 points are found, design with best channel aspect ratio must be found

    def colebrook(self, D, Re, tol=1e-6, max_iter=100):
        lambda_f = 0.02
        for _ in range(max_iter):
            rhs = -0.869 * np.log((self.epsilon / (3.7 * D)) + (2.523 / (Re * np.sqrt(lambda_f))))
            lambda_new = 1 / (rhs ** 2)
            if abs(lambda_new - lambda_f) < tol:
                return lambda_new
            lambda_f = lambda_new
        raise RuntimeError("Colebrook did not converge")
        
    def compute_pressure_drop(self, lambda_f, D, U):
        return lambda_f * (self.chanel_length / D) * (self.density * U**2 / 2)

    def find_minimum_diameter(self, n, dp_max=1e5, D_start=0.1e-3, D_stop=100e-3, D_step=0.01e-3):
        for D in np.arange(D_start, D_stop, D_step):
            A = np.pi * (D**2) / 4
            U = (self.injector.fuel_massflowrate / n) / (self.density * A)
            # if U > 1:
            #     print(U, D*1000)
            Re = (self.density * U * D) / self.dynamic_viscosity
            if Re < 4000:
                continue  # skip laminar regions
            try:
                lambda_f = self.colebrook(D, Re)
            except RuntimeError:
                continue
            dp = self.compute_pressure_drop(lambda_f, D, U)
            if dp <= dp_max:
                return D
        # raise ValueError(f"No suitable diameter found in given range: D_start={D_start}, D_stop={D_stop}, D_step={D_step}, nchannesl={n}")
        return None

    def print_results(self):
        print("---------------------------Optimisation results:----------------------------")
        print("Metal thickness domain (m):", self.metal_thickness)
        print("Number of channels domain:", self.number_of_channels)
        print("Channel area domain (m^2):", self.povrsina_kanala)

        # Visualize the full 3D domain (all combinations)
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(
                self.domain_points[:, 0],  # Metal thickness
                self.domain_points[:, 1],  # Number of channels
                self.domain_points[:, 2],  # Channel area
                c='b', marker='o', s=2
            )
            ax.set_xlabel('Metal Thickness [m]')
            ax.set_ylabel('Number of Channels')
            ax.set_zlabel('Channel Area [m^2]')
            ax.set_title('3D Optimisation Domain (Full Grid)')
            plt.show()
        except ImportError:
            print("matplotlib is not installed, skipping 3D plot.")

        #if self.calculate_minimum_area() > self.povrsinaKanala(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, self.fi_inlet_rad_ulaz):
        #    print("channel area too small")



