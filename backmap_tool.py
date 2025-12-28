if add_nonbonded_force:
    if nonbonded_type == "gaussian":
        # Stage 1: Gaussian "Bulldozer" Force
        # Used for initial untangling. Ignores chemistry, only geometric repulsion.
        # Formula: E = Height * exp( - (r / Width)^2 )
        self.nb_force = mm.CustomNonbondedForce(
            "ga_k * ga_h * exp(-(r/ga_w)^2);"
        )
        # Global parameters (can be tuned in Context)
        self.nb_force.addGlobalParameter("ga_k", 1.0)    # ON/OFF switch
        self.nb_force.addGlobalParameter("ga_h", 500.0)  # Height (kJ/mol)
        self.nb_force.addGlobalParameter("ga_w", 0.1)    # Width (nm)
        
        # Dummy parameters to match the addParticle signature in _set_nonbonded_parameters
        # _set_nonbonded_parameters calls addParticle([atomType, q])
        self.nb_force.addPerParticleParameter("dummy_type") 
        self.nb_force.addPerParticleParameter("dummy_q")
        
        self.nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        self.nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
        sys.addForce(self.nb_force)
        
        # Note: In Gaussian mode, we do NOT add exception forces (es/lj_except_force)
        # because we want pure geometric repulsion to untangle knots.
        
    elif nonbonded_type == "softcore":
        # Stage 2: Gapsys Soft-Core Force (C6/C12 Version)
        # Uses Tabulated Functions for C6/C12 but modifies the expression for softness.
        energy_expression = (
            "step(rcut - r) * (LJ_soft - corrLJ + ES_soft);"
            
            # Soft-Core LJ (modified Gapsys form using C12/C6 directly)
            # V = C12 / (reff^6)^2 - C6 / reff^6
            "LJ_soft = (C12(type1, type2) / reff_sq) - (C6(type1, type2) / reff);"
            "reff_sq = reff * reff;"
            
            # Softened distance: reff^6 = r^6 + alpha * (C12/C6)
            # Use 'select' to handle C6=0 (pure repulsion) cases safely
            "reff = r^6 + soft_alpha * sigma6_proxy;"
            "sigma6_proxy = select(step(1e-8 - abs(C6(type1, type2))), 1.0, C12(type1, type2) / C6(type1, type2));"
            
            # Soft-Core Electrostatics
            "ES_soft = f * q1 * q2 / epsilon_r * (inv_r_soft + krf * r^2 - crf);"
            "inv_r_soft = 1 / sqrt(r^2 + soft_alpha_coul);"
            
            # Long Range Correction (Switch)
            "corrLJ = step(r - rswitch) * (C12(type1, type2)/r^12 - C6(type1, type2)/r^6) * (10 * x^3 - 15 * x^4 + 6 * x^5);"
            "x = (r - rswitch) / (rcut - rswitch);"
            
            # Constants
            "rswitch = 0.9;"
            f"rcut = {nonbonded_cutoff.value_in_unit(unit.nanometers)};"
            "krf = 1 / (2 * rcut^3);"
            "crf = 1 / rcut + krf * rcut^2;"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
        )
        
        self.nb_force = mm.CustomNonbondedForce(energy_expression)
        
        # Parameters
        self.nb_force.addGlobalParameter("soft_alpha", 0.5)      # VDW Softness
        self.nb_force.addGlobalParameter("soft_alpha_coul", 0.1) # Electrostatic Softness
        
        self.nb_force.addPerParticleParameter("type")
        self.nb_force.addPerParticleParameter("q")
        
        self.nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        self.nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
        sys.addForce(self.nb_force)
        
    else: # "standard" (Original code logic)
        self.nb_force = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corrLJ + ES);"
            "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
            "corrLJ = step(r-rswitch) * (C12(type1, type2) / r^12 - C6(type1, type2) / r^6) * (10 *x^3 - 15 * x^4 + 6 * x^5);"
            "x = (r - rswitch)/(rcut - rswitch);"
            "rswitch = 0.9;"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        #self.nb_force.addComputedValue('k', 'step(0.5 - soft) * 1 + step(soft - 0.5) * (1 + cos(3.14159265359 * r / rcut)) ')
        self.nb_force.addGlobalParameter("soft", 1)
        self.nb_force.addPerParticleParameter("type")
        self.nb_force.addPerParticleParameter("q")
        self.nb_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        self.nb_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))
        sys.addForce(self.nb_force)