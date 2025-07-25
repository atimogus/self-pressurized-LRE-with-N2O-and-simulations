import numpy as np

class regenNozzle:
    def __init__(self, combustion_chamber_diameter, vanjski_promjer_ulaz, promjer_grla_mlaznice, tlak_komore, tlak_fluida, cvrstoca_materijala, faktor_sigurnosti, 
                 d_hidraulicke_cijevi, maseni_protok_IPA, gustoca_IPA, alu_thickness, n, k=1):
        # Input parameters
        self.vanjski_promjer_ulaz = vanjski_promjer_ulaz  # mm
        self.promjer_grla_mlaznice = promjer_grla_mlaznice  # mm
        self.tlak_komore = tlak_komore  # bar
        self.tlak_fluida = tlak_fluida  # bar
        self.cvrstoca_materijala = cvrstoca_materijala  # MPa
        self.faktor_sigurnosti = faktor_sigurnosti
        self.d_hidraulicke_cijevi = d_hidraulicke_cijevi  # mm
        self.maseni_protok_IPA = maseni_protok_IPA  # kg/s
        self.gustoca_IPA = gustoca_IPA  # kg/m^3
        self.alu_thickness = alu_thickness  # mm
        self.n = n  # number of channels
        self.k = k  # expansion ratio
        self.combustion_chamber_diameter = combustion_chamber_diameter  # mm

        # Derived parameters
        self.Auk = np.pi * (self.d_hidraulicke_cijevi / 1000)**2 / 4  # m^2

        self.fi_alu_rad, self.fi_alu_stupnjevi = self.alu_stjenka_geometry()

        """ Calculate dimensioons for cooling channels on the throat of the nozzle """
        self.inner_wall_thickness_throat = self.calculate_wall_thickness(self.tlak_fluida, self.promjer_grla_mlaznice/1000)
        self.d_cooling_channels_throat = self.promjer_grla_mlaznice / 1000 + 2 * self.inner_wall_thickness_throat
        self.D_cooling_channels_throat = self.calculate_external_diameter(self.d_cooling_channels_throat, self.fi_alu_rad)
        self.outer_wall_thickness_throat = self.calculate_wall_thickness(self.tlak_fluida , self.D_cooling_channels_throat)
        self.D_throat = self.D_cooling_channels_throat + 2 * self.outer_wall_thickness_throat

        """ Calculate dimensioons for cooling channels on the inlet of the nozzle """
        self.inner_wall_thickness_inlet = self.calculate_wall_thickness(self.tlak_komore, self.combustion_chamber_diameter)
        self.d_cooling_channels_inlet = self.combustion_chamber_diameter + 2 * self.inner_wall_thickness_inlet
        self.D_cooling_channels_inlet = self.calculate_external_diameter(self.d_cooling_channels_inlet, self.fi_alu_rad)
        self.outer_wall_thickness_inlet = self.calculate_wall_thickness(self.tlak_fluida , self.D_cooling_channels_inlet)
        self.D_inlet = self.D_cooling_channels_inlet + 2 * self.outer_wall_thickness_inlet
        self.channel_height = self.channel_height()

        """ Calculate all angles for the inlet and throat of the nozzle """
        self.fi_inlet_rad_ulaz, self.fi_inlet_stupnjevi_ulaz = self.calculate_fi_ulaz(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet)
        self.fi_inlet_rad_izlaz, self.fi_inlet_stupnjevi_izlaz = self.calculate_fi_izlaz(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet)
        self.fi_grlo_rad_ulaz, self.fi_grlo_stupnjevi_ulaz = self.calculate_fi_ulaz(self.D_cooling_channels_throat, self.d_cooling_channels_throat)
        self.fi_grlo_rad_izlaz, self.fi_grlo_stupnjevi_izlaz = self.calculate_fi_izlaz(self.D_cooling_channels_throat, self.d_cooling_channels_throat)
        self.fi_revolve = self.fi_inlet_stupnjevi_ulaz + self.fi_inlet_stupnjevi_izlaz + 2 * self.fi_alu_stupnjevi

        self.minimal_channel_lenght = self.minimal_channel_lenght()
        self.brzina_ucijevi = self.brzina_protoka_uCijevi()
        self.brzina_uKanalu_grlo = self.brzina_protoka_uKanalu(self.D_cooling_channels_throat, self.d_cooling_channels_throat, self.fi_grlo_rad_ulaz)
        self.brzina_uKanalu_inlet = self.brzina_protoka_uKanalu(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, self.fi_inlet_rad_ulaz)
        if self.brzina_uKanalu_grlo - self.brzina_uKanalu_inlet > 0.000001:
            raise ValueError("Calculated velocities mismatch detected.")

        # test dizajna kanala
        ukupna_povrsina_test = (self.n * ( self.povrsinaKanala(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, self.fi_inlet_rad_ulaz) +
                                                                   self.povrsinaKanala(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, self.fi_inlet_rad_izlaz)) + 
                                                                   2 * self.n * self.povrsinaKanala(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, self.fi_alu_rad))
        ukupna_povrsina_inlet = self.povrsinaKanala(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, 2 * np.pi)
        if ukupna_povrsina_test - ukupna_povrsina_inlet > 0.000001:
            raise ValueError("Calculated surface area mismatch detected.")

    def calculate_wall_thickness(self, tlak, promjer): #tlak mora biti u paskalima, promjer u metrima
        """Calculate wall thickness"""
        t1 = max(np.ceil((self.faktor_sigurnosti * tlak * promjer) / 
                 (2 * self.cvrstoca_materijala) * 1e5) / 1e5, 1.2 * 1e-3)  
        t2 = max(np.ceil((self.faktor_sigurnosti * tlak * promjer) / 
                 (4 * self.cvrstoca_materijala) * 1e5) / 1e5, 1.2 * 1e-3) 
        # print(f"t1: {t1:.6f} m")
        # print(f"t2: {t2:.6f} m")
        if t1 > t2:
            return t1
        else:
            return t2
    
    def calculate_external_diameter(self, d, fi_rad):
        """Calculate external diameter wit internal diameter."""
        D = np.sqrt(4 / np.pi * self.Auk * (1 + self.k) / (1 - self.n * fi_rad / np.pi)  + d**2)  
        return D
    
    def calculate_internal_diameter(self, D, fi_rad):
        """Calculate internal diameter with external diameter."""
        d = np.sqrt(D**2 - (4 / np.pi * self.Auk * (1 + self.k))/(1 - self.n * fi_rad/(np.pi)))
        return d

    def alu_stjenka_geometry(self):
        """Calculate the geometry of the aluminum wall that separates inlet and outlet."""
        wall_thickness = self.calculate_wall_thickness(self.tlak_fluida - self.tlak_komore, self.promjer_grla_mlaznice/1000)
        d_rashladnih_kanala = self.promjer_grla_mlaznice / 1000 + 2 * wall_thickness
        fi_alu_rad = 2 * (self.alu_thickness) / (d_rashladnih_kanala)
        fi_alu_stupnjevi = np.degrees(fi_alu_rad)
        return fi_alu_rad, fi_alu_stupnjevi


    def calculate_fi_ulaz(self, D, d):
        """Calculate fi_rad and fi_deg for a given geometry.    """
        # Calculate fi_rad and fi_deg for the inlet
        fi_ulaz_rad = self.Auk * (8 / self.n) * (1 / (D**2 - d**2))
        fi_ulaz_deg = np.degrees(fi_ulaz_rad)

        return fi_ulaz_rad,fi_ulaz_deg
    
    def calculate_fi_izlaz(self, D, d):
        """Calculate fi_rad and fi_deg for a given geometry.    """

        # Calculate fi_rad and fi_deg for the outlet
        fi_izlaz_rad = self.k * self.Auk * (8 / self.n) * (1 / (D**2 - d**2))
        fi_izlaz_deg = np.degrees(fi_izlaz_rad)

        return fi_izlaz_rad, fi_izlaz_deg

    def brzina_protoka_uCijevi(self):
        """Calculate the velocity of IPA."""
        brzina_toka_IPA = self.maseni_protok_IPA / (self.gustoca_IPA * self.Auk)
        return brzina_toka_IPA

    def brzina_protoka_uKanalu(self , D , d , fi_rad):
        """Calculate the velocity of IPA."""
        povrsinaKanala = (D**2 -d**2) * fi_rad/8
        brzina_toka_IPA = (self.maseni_protok_IPA/self.n) / (self.gustoca_IPA * povrsinaKanala)
        return brzina_toka_IPA

    def povrsinaKanala(self, D, d, fi_rad):
        """Calculate the area of the channel."""
        return (D**2 -d**2) * fi_rad/8
    
    def channel_height(self):
        visina_prolaz_izmedju_kanala = 3 * (self.Auk / 
                                            (self.n * 
                                             (self.D_cooling_channels_throat - self.d_cooling_channels_throat)))
        return visina_prolaz_izmedju_kanala
    
    def minimal_channel_lenght(self):
        """Calculate the minimal channel length."""
        # opseg_d_cooling_channels_throat *  self.fi_grlo_rad_ulaz
        minimal_length_thickness = np.pi * self.d_cooling_channels_throat * self.fi_grlo_rad_ulaz
        minimal_lenght_radius  = self.D_cooling_channels_inlet/2 - self.d_cooling_channels_inlet/2

        if minimal_length_thickness < minimal_lenght_radius:
            minimal_length = minimal_length_thickness
        else:
            minimal_length = minimal_lenght_radius

        return minimal_length


    def printDesign(self):
        """Print the calculated values."""
        print("\n ------------------------REGEN NOZZLE GEOMETRY-------------------------------------------")
        
        #nozzle throat
        print(f"promjer_grla_mlaznice: {self.promjer_grla_mlaznice:.2f} mm")
        print(f"d_cooling_channels_throat: {self.d_cooling_channels_throat * 1000:.2f} mm")
        print(f"D_cooling_channels_throat: {self.D_cooling_channels_throat * 1000:.2f} mm")
        print(f"D_throat: {self.D_throat * 1000:.2f} mm\n")

        #combustion chamber
        print(f"Combustion Chamber Diameter: {self.combustion_chamber_diameter * 1000:.2f} mm")
        print(f"inner_wall_thickness_inlet: {self.inner_wall_thickness_inlet * 1000:.2f} mm")
        print(f"d_cooling_channels_inlet: {self.d_cooling_channels_inlet * 1000:.2f} mm")
        print(f"D_cooling_channels_inlet: {self.D_cooling_channels_inlet * 1000:.2f} mm")
        print(f"D_inlet: {self.D_inlet * 1000:.2f} mm\n")

        print(f"fi_alu_rad: {self.fi_alu_rad:.6f} rad, fi_alu_stupnjevi: {self.fi_alu_stupnjevi:.6f}°")
        print(f"fi_inlet_rad_ulaz: {self.fi_inlet_rad_ulaz:.6f} rad = {self.fi_inlet_stupnjevi_ulaz:.6f}°")
        print(f"fi_inlet_rad_izlaz: {self.fi_inlet_rad_izlaz:.6f} rad = {self.fi_inlet_stupnjevi_izlaz:.6f}°\n")
        print(f"fi_grlo_rad_ulaz: {self.fi_grlo_rad_ulaz:.6f} rad = {self.fi_grlo_stupnjevi_ulaz:.6f}°")
        print(f"fi_grlo_rad_izlaz: self.{self.fi_grlo_rad_izlaz:.6f} rad = {self.fi_grlo_stupnjevi_izlaz:.6f}°\n")
        print(f"kut za koj je potrebno napraviti revolve: {self.fi_revolve:.6f}°")

        print(f"Visina prolaza izmedju kanala: {self.channel_height * 1000:.2f} mm")

        print(f"brzina_ucijevi: {self.brzina_ucijevi:.6f} m/s")
        print(f"ukupna_povrsina_inlet-ukupna za sve unutarnje kanale: {self.n * ( self.povrsinaKanala(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, self.fi_inlet_rad_ulaz)):.8f} m^2")
        print(f"ukupna_povrsina_inlet za jedan kanal: {self.povrsinaKanala(self.D_cooling_channels_inlet, self.d_cooling_channels_inlet, self.fi_inlet_rad_ulaz):.8f} m^2")



# # Example usage
# if __name__ == "__main__":
#     nozzle = NozzleDesign(
#         vanjski_promjer_ulaz=160,
#         promjer_grla_mlaznice=2 * 16.7,
#         tlak_komore=30e5,
#         tlak_fluida=70e5,
#         cvrstoca_materijala=230e6,
#         faktor_sigurnosti=1.5,
#         d_hidraulicke_cijevi=25,
#         maseni_protok_IPA=0.4,
#         gustoca_IPA=785,
#         alu_thickness=2e-3,
#         n=12,
#         k=1
#     )

#     # Calculate the geometry of the aluminum wall
#     fi_alu_rad, fi_alu_stupnjevi = nozzle.alu_stjenka_geometry()
#     print(f"fi_alu_rad: {fi_alu_rad:.6f} rad")
#     print(f"fi_alu_stupnjevi: {fi_alu_stupnjevi:.6f}°\n")

#     # Calculate the geometry of the channels on the throat
#     inner_wall_thickness_throat = nozzle.calculate_wall_thickness(nozzle.tlak_fluida - nozzle.tlak_komore, nozzle.promjer_grla_mlaznice/1000)
#     d_cooling_channels_throat = nozzle.promjer_grla_mlaznice / 1000 + 2 * inner_wall_thickness_throat
#     D_cooling_channels_throat = nozzle.calculate_external_diameter(d_cooling_channels_throat, fi_alu_rad)
#     outer_wall_thickness_throat = nozzle.calculate_wall_thickness(nozzle.tlak_fluida , D_cooling_channels_throat)
#     D_throat = D_cooling_channels_throat + 2 * outer_wall_thickness_throat

#     print(f"D_cooling_channels_throat: {D_cooling_channels_throat * 1000:.2f} mm\n")
#     print(f"D_throat: {D_throat * 1000:.2f} mm\n")

#     # Calculate the geometry of the channels on the nozzle inlet
