# -*- coding: utf-8 -*-
""" 
Default parameters for model setup and payoff metrics
    1. Changes here is cause issues with loading previous simulations, compatibility not assured
    2. These values may be modified by the solver as needed to obtain valid solutions
"""
#initializations
import numpy as np
if __package__ is None or __package__ == '':
    import properties
    import stress
    import stats
    import csv
    from units import *
else:
    from . import properties
    from . import stress
    from . import csv
    from . import stats
    from .units import *
water = properties.water()
deg=deg
yr=yr
cP=cP
darcy=darcy
mD=mD
g=g

#**********************************************************************************************************
### bulk transfer of scalar class attributes
#**********************************************************************************************************
def replace_attributes(source,target):
    #match fields from setup to cashflow
    for item in vars(source):
        if item in vars(target):
            setattr(target,item,getattr(source,item))

def transfer_attributes(source,target):
    #store attributes with match fields
    for item in vars(source):
        if np.isscalar(getattr(source,item)):
            setattr(target,item,getattr(source,item))

#**********************************************************************************************************
### independent parameters
#**********************************************************************************************************
class setup:
    def __init__(self,run=True):
        #identifier
        self.PIN_PIN_PIN = np.random.randint(10000000,99999999)
        
        #GeoDT solver parameters
        self.Strategy_Iterations_units = 10
        self.Strategy_DiscreteFractures_path = '' #os.getcwd()+'\\dfn.csv'
        self.Strategy_DiscreteWells_path = '' #os.getcwd()+'\\well.csv'
        self.Strategy_Solver_type = 'gg' #'Gringarten' 'GeoDT' 'gg' 'gt'
        self.Strategy_Design_type = 'O&G' #'ALL' 'EGS' 'AGS' 'CGS' 'O&G' 'DCM' 'FGS'
        self.Strategy_Target_type = 'Deep' #'Deep','Temp'
        self.Strategy_Stim_type = 'Radi' #'Volu','Radi'
        self.Strategy_Cluster_type = 'Num' #'Var','Count','Number','Spa','Spacing'
        self.Strategy_SecondOrderNoise_ratio = 0.0 #0.0 = no second order randomness, 1.0 = full, 2.0 = amplified #TODO
        self.Well_TargetTemp_K = 273.15+325 #K
        self.Well_TargetDepth_m = 6000.0 #m
        self.Well_ProducerCount_wells = 2 #wells
        self.Well_Pattern_count = 1 #pattern #TODO: needs to be implemented
        self.Well_Spacing_m = 175.0 #m
        self.Well_DeviatedLength_m = 3000.0 #m
        self.Well_ProducerProportion_ratio = 0.8 #m/m
        self.Well_Azimuth_rad = 90.0*deg #rad
        self.Well_Dip_rad = 30.0*deg #rad
        self.Well_RotationPhase_rad = 0.0*deg #rad
        self.Well_RotationToe_rad = 0.0*deg #rad
        self.Well_RotationSkew_rad = 0.0*deg #rad
        self.Well_ProductionDiameter_m = 0.0254*8.625 #m #TODO: Implement complete well plan
        self.Well_WallThickness_m = 0.0254*0.75 #m #TODO: Implement complete well plan
        self.Well_AnnularThickness_m = 0.0254*1.25 #m #TODO: Implement complete well plan
        self.Well_SegmentCount_units = 10.0 #units, number of breaks in the well for AGS, first is conductor #TODO: Implement complete well plan
        self.Well_Roughness_metric = 80.0
        self.Circulation_TargetRate_m3ps = 0.10 #m3/s, flow rate per well
        self.Stimulation_InjectionRate_m3ps = 0.40 #m3/s, flow rate per well
        self.Stimulation_TargetVolume_m3pfrac = 50000.0 #m3 #TODO: update volume calculations to m3 per fracture
        self.Stimulation_TargetRadius_m = 196.0 #m3
        self.Stimulation_Leakoff_ratio = 0.1 #3.55 #m3/m3, volume leaked vs volume injected (e.g., 0.1 = 10% loss)
        self.Circulation_Leakoff_ratio = 0.1 #3.55 #m3/m3, volume leaked vs volume injected (e.g., 0.1 = 10% loss)
        self.Stimulation_CriticalPressure_Pa = 2e6 #Pa
        self.Stimulation_Steps_count = 6 #steps
        self.Domain_Size_m = 1000.0 #m
        self.Rock_ThermalGradient_Kpm = 0.042 #K/m
        self.Rock_Density_kgpm3 = 2700.0 # kg/m3
        self.Rock_ThermalConductivity_WpmK = 2.5 # W/m-K
        self.Rock_HeatCapacity_kJpm3K = 2063.0 # kJ/m3-K
        self.Rock_Biot_ratio = 0.0 #Pa/Pa
        self.Air_Temperature_K = 273.15 # K
        self.Air_Pressure_Pa = 0.101*1e6 #Pa
        self.Rock_YoungsModulus_Pa = 50.0*1e9 #Pa
        self.Rock_PoissonRatio_ratio = 0.3 #
        self.Stress_ShminToSv_ratio = 0.5 #Pa/Pa
        self.Stress_ShmaxToShmin_ratio = 1.5 #Pa/Pa
        self.Stress_ShminAzimuth_rad = 115.0*deg #rad
        self.Stress_ShminDip_rad = 15.0*deg #rad
        self.Stress_Variance_rad = 0.5*deg #rad
        self.Shear_FrictionCoefficient_ratio = 0.6 #T:N #TODO
        self.Shear_Cohesion_Pa = 2e6 #Pa #TODO
        self.Shear_SlipLength_ratio = 0.010 #m/m #TODO
        self.Shear_DilationSlip_ratio = 0.100 #m/m #TODO
        self.Shear_HydraulicDilation_ratio = 0.6 #m/m #TODO
        self.Tension_FrictionCoefficient_ratio = 0.56 #rad #TODO
        self.Tension_Cohesion_Pa = 0.2*1e6 #Pa #TODO
        self.Tension_FractureCompressibility_1pPa = 2.9e-8 #1/Pa #TODO
        self.Tension_Toughness_Pasqrtm = 1.5*1e6 #Pa-m**0.5
        self.Stimulation_ProppantDensity_kgpm3 = 1538.0 #kg/m3 #TODO: solid density (no porosity)
        self.Stimulation_ProppantCompressibility_1pPa = 2.9e-8 #1/Pa #TODO
        self.JointSets_InitialDilation_m = 0.00010 #m #TODO
        self.Domain_BoundaryAperture_m = 0.003 #m
        self.Fracture_Tortuosity_ratio = 0.90 #m/m #TODO
        self.Cement_ThermalConductivity_WpmK = 2.0 # W/m-K
        self.Cement_HeatCapacity_kJpm3K = 2000.0 # kJ/m3-K
        self.Power_GeneralEfficiency_ratio = 0.85 # kWe/kWt
        self.Power_LifeSpan_s = 30*yr #years
        self.Domain_TimeSteps_steps = 30 #steps
        self.Circulation_Backpressure_Pa = 10.0*1e6 #Pa
        self.Circulation_InjectionTemperature_K = 60.0+273.15 #K
        self.Power_ConvectionCoefficient_kWpm2K = 3.0 #kW/m2-K
        self.Rock_Permeability_m2 = 0.1*mD #m2
        self.Stimulation_ProppantPermeability_m2 = 100.0*darcy #m2 #TODO
        self.JointSets_Set1_fractures = 16 #fractures
        self.JointSets_Set2_fractures = 8 #fractures
        self.JointSets_Set3_fractures = 4 #fractures
        self.JointSets_Variance_rad = 7.0*deg #rad
        self.JointSets_DiameterMin_m = 10.0 #m
        self.JointSets_DiameterMax_m = 1000.0 #m
        #self.Stimulation_Clusters_clusters = 81 #clusters
        self.Stimulation_TargetClusters_count = 81 #clusters
        self.Stimulation_TargetSpacing_m = 3000.0/(81+1) #frac spacing
        self.Stimulation_PerfDiameter_m = 0.013 #m
        self.Stimulation_PerfPerCluster_perfs = 6 #holes
        self.Stimulation_ProppantConcentration_m3pm3 = 0.045 #m3/m3
        self.Circulation_AllowableOverpressure_Pa = 0.99 #pinj/s3 #Pa

        #Drilling Cost Model
        self.Drilling_Utilization_ratio = 1.00 #active use vs idle time for drilling rig
        self.Drilling_TimeMobilize_h = 336.0 #hr
        self.Drilling_CostMobilize_USD = 2000000.0 #5130000.0 #USD
        self.Drilling_TimeRigWalk_h = 48.0 #hr
        self.Drilling_CostRigWalk_USD = 250000.0 #USD
        self.Drilling_TimeCut_hpm2 = 6.6e-9 #h/m2
        self.Drilling_TimeCut_hpm = 0.020 #h/m
        self.Drilling_CostCut_USDph = 4000.0 #5082.1 #$/hr
        self.Drilling_CostCut_USDpm4 = 1.0e-5 #2.0e-5 #$/m4
        self.Drilling_CostCut_USDpm3 = 1500.0 #1000.0 #$/m3
        self.Drilling_TimeCasing_hpm = 0.01 #h/m
        self.Drilling_CostCasing_USDph = 3500.0 #4001.3 #$/hr
        self.Drilling_CostCasing_USDpm4 = 0.5 #1.2 #$/m4
        self.Drilling_CostCasing_USDpm3 = 1500.0 #2000.0 #$/m3
        self.Drilling_Learning_ratio = 0.30 #0.20 #learning rate
        
        #Stimulation Cost Model
        self.Fracking_Fixed_USD = 500000.0 #1233149.0 #$
        self.Fracking_Depth_USD = 1000.0 #638.0 #$/m
        self.Fracking_Hourly_USDph = 5000.0 #5973.3 #$/hr
        self.Fracking_Sand_USDpkg = 0.662 #$/kg #ISP
        self.Fracking_Fuel_USDpkWh = 0.15 #0.451 #$/kWh
        self.Fracking_Water_USDpm3 = 22.2 #$/m3
        self.Fracking_PressureFactor_ratio = 0.8
        self.Fracking_PressureFactor_scale = 1e-8 #Pa
        self.Fracking_EquipmentWear_USDpm3 = 9.85 #$/m3
        self.Fracking_TimeFactor_hph = 1.95
        
        #Economic Parameters
        self.Cost_Planning_USD = -3.5e6 #USD, survey & planning costs
        self.Cost_PlanningTime_yr = 0.6 #yr, site prep and well planning
        self.Cost_ConstructionTime_yr = 2.0 #yr, minimum surface construction time
        self.Cost_Gridconnect_USD = -4.0e6 #USD, grid connection costs
        self.Cost_EquipmentCost_USDpkW = 2500.0 #%/kWe
        self.Cost_Grants_USD = 0.0e6 #USD, government grants
        self.Cost_Operations_ratio = 0.02 #ratio, opex vs capex
        self.Cost_Seismic_USDpMw = 2e-4 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re)
        self.Cost_Seismic_exp = 5.0 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re)
        self.Sales_ElectricWholesale_USDpkWh = 54.645e-3 #USD, market rate sales for electricity
        self.Sales_ElectricContract_USDpkWh = 165.590e-3 #USD, contract sales rate for electricity
        self.Sales_ElectricPurchase_USDpkWh = 291.300e-3 #USD, typical purchase price for customers
        self.Sales_HeatWholesale_USDpkWh = 10.000e-3 #USD, market rate sales for heat
        self.Sales_HeatContract_USDpkWh = 10.000e-3 #USD, contract rate sales for heat
        self.Sales_HeatDemand_kWhpyr = 100578.8571 #kWh per yr, demand for heat delivery
        self.Sales_CapacityFactor_ratio = 0.95 #ratio, capacity factor for plant uptime
        self.Sales_ContractTime_yr = 15 #yr, length of contract for heat and electricity
        self.Equity_CostOfCapital_ratio = 0.12 #ratio, time value of money versus alternative investments
        self.Equity_Contribution_ratio = 0.40 #ratio, equity contribution assuming debt service reserve account funded by debt 
        self.Equity_CompanyShare_ratio = 0.02 #ratio, share of equity held by primary company
        self.Equity_SeriesAShare_ratio = 0.08 #ratio, share of equity held by series A investors
        self.Equity_SeriesBShare_ratio = 0.90 #ratio, share of equity held by series B investors
        self.Equity_BrokerFee_ratio = 0.050 #ratio, fee to setup project funds
        self.Debt_InterestRate_ratio = 0.080 #ratio, interest rate on debt and loans 
        self.Debt_ArrangementFee_ratio = 0.015 #ratio, fee to setup debt accounts
        self.Debt_ServiceAccount_ratio = 0.50 #principal and interest fraction of a year
        self.Debt_Tenor_yr = 30.0 #years, loan duration until maturity
        self.Tax_Inflation_ratio = 0.02 #inflation rate
        self.Tax_Corporate_ratio = 0.40 #ratio, tax rate on profits after losses 
        self.Tax_ProductionCredit_USDpkWh = 0.0 #0.015 #USD, tax credit for energy production
        self.Tax_CapitalCredit_ratio = 0.0  #0.40 #USD, tax credit for capital construction
        self.Tax_DirectPay_ratio = 0.0 #1.0 #ratio, ratio of tax credits that are elligble for direct pay (e.g., Inflation Reduction Act)
        self.Tax_DepreciationRate_ratio = 0.20
        self.CO2_Electric_kgpkWh = 0.22499 #kgCO2/kWhe, CO2 emitted from 'grid' power
        self.CO2_Heat_kgpkWh = 0.18293 #kgCO2/kWht, CO2 emitted from natural gas for heat 
        self.CO2_Geothermal_kgpkWh = 0.00154 #kgCO2/kWht, CO2 emitted from natural gas for heat
        
        #*******************************************************************************************************************************************************
        """
        #relics for removal #!!! remove these before release
        self.Domain_PressureIncrement_Pa = 0.1*1e6 #Pa
        self.Economics_ElectricitySales_USDpkWh = 0.1372 #$/kWh - customer electricity retail price 
        self.Economics_DrillingCost_USDpm = 2763.06 #$/m - Lowry et al, 2017 large diameter well baseline 
        self.Economics_AccessCost_USD = 590e3 #$ Lowry et al, 2017 large diameter well baseline 
        self.Economics_EquipmentCost_USDpkW = 2025.65 #$/kWe simplified from GETEM model 
        self.Economics_ExplorationCost_USDpm = 2683.41 #$/m simplified from GETEM model 
        self.Economics_Maintenance_USDpkWh = 0.03648 #$/kWh simplified from GETEM model 
        self.Economics_Seismic_USDpMw = 2e-4 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re) 
        self.Economics_Seismic_exp = 5.0 #$/Mw for $300M Mw 5.5 quake Pohang (Westaway, 2021) & $17.2B Mw 6.3 quake Christchurch (Swiss Re) 
        self.Well_HydraulicRadius_m = 0.0254*5.0 #m 
        self.Well_CasingRadius_m = 0.0254*5.5 # m 
        self.Well_BoreRadius_m = 0.0254*6.5 # m 
        self.StressMin_Coefficient_ratio = 0.5 #Pa/Pa 
        self.StressMax_Coefficient_ratio = 0.75 #Pa/Pa 
        self.StressMin_Azimuth_rad = 90.0*deg #rad 
        self.StressMin_Dip_rad = 0.0*deg #rad 
        self.StressMin_Variance_rad = 0.5*deg #rad
        self.Fracture_CompressibilityMin_1pPa = 2.0e-9 #1/Pa
        self.Fracture_CompressibilityNom_1pPa = 2.9e-8 #1/Pa
        self.Fracture_CompressibilityMax_1pPa = 10.0e-8 #1/Pa
        self.Proppant_CompressibilityMin_1pPa = 2.0e-9 #1/Pa
        self.Proppant_CompressibilityNom_1pPa = 2.9e-8 #1/Pa
        self.Proppant_CompressibilityMax_1pPa = 10.0e-8 #1/Pa
        self.Proppant_InitialPermeabilityMin_m2 = 10.0*darcy #m2
        self.Proppant_InitialPermeabilityNom_m2 = 100.0*darcy #m2
        self.Proppant_InitialPermeabilityMax_m2 = 300.0*darcy #m2
        self.Fracture_TortuosityMin_ratio = 0.80 #m/m
        self.Fracture_TortuosityNom_ratio = 0.90 #m/m
        self.Fracture_TortuosityMax_ratio = 1.00 #m/m
        self.Fluid_NominalDensity_kgpm3 = 965.0 #kg/m3
        self.Fluid_Viscosity_Pas = 0.2*cP #Pa-s
        self.Fracture_SlipLengthMin_ratio = 10.0**-3.0 #m/m
        self.Fracture_SlipLengthNom_ratio = 10.0**-2.0 #m/m
        self.Fracture_SlipLengthMax_ratio = 10.0**-1.2 #m/m
        self.Fracture_DilationSlipMin_ratio = 0.000 #m/m
        self.Fracture_DilationSlipNom_ratio = 0.200 #m/m
        self.Fracture_DilationSlipMax_ratio = 0.800 #m/m
        self.Fracture_HydraulicDilationMin_ratio = 0.0 #m/m
        self.Fracture_HydraulicDilationNom_ratio = 0.6 #m/m
        self.Fracture_HydraulicDilationMax_ratio = 2.0 #m/m
        self.Fracture_InitialHydraulicApertureMin_m = 0.00005 #m
        self.Fracture_InitialHydraulicApertureNom_m = 0.00010 #m 
        self.Fracture_InitialHydraulicApertureMax_m = 0.00020 #m 
        self.Shear_FrictionAngleMin_rad = 20.0*deg #rad
        self.Shear_FrictionAngleNom_rad = 35.0*deg #rad
        self.Shear_FrictionAngleMax_rad = 45.0*deg #rad
        self.Shear_CohesionMin_Pa = 1.0*1e6 #Pa
        self.Shear_CohesionNom_Pa = 7.0*1e6 #Pa
        self.Shear_CohesionMax_Pa = 15.0*1e6 #Pa
        self.Domain_StimulationLimit_steps = 1
        self.Strategy_SourceDirectory_path = os.getcwd()
        self.Strategy_TargetDirectory_path = os.getcwd()
        self.Strategy_Completion_type = 'DC' #type of completion; 'DC' = direct contact 'FE' = fracture exchange
        #*******************************************************************************************************************************************************
        """
        
        if run: self.re_init()
        
#**********************************************************************************************************
### dependent parameters (e.g., stress, number of injectors, temperature)
#**********************************************************************************************************
    def re_init(self):
        #get new pin
        self.PIN_PIN_PIN = np.random.randint(10000000,99999999)
        
        '''
        self.Stimulation_TargetClusters_count = 81 #clusters
        self.Stimulation_TargetSpacing_m = 3000.0/(81+1) #frac spacing
        '''
        
        #fracture spacing
        if self.Strategy_Cluster_type.upper() in ['NUM','FIX','COUNT','NUMBER']:
            if self.Stimulation_TargetClusters_count > 1:
                self.Stimulation_TargetSpacing_m = self.Well_DeviatedLength_m*self.Well_ProducerProportion_ratio/(self.Stimulation_TargetClusters_count+1)
            else:
                self.Stimulation_TargetSpacing_m = self.Domain_Size_m
        else: #in ['VAR','SPA','SPACING']
            self.Stimulation_TargetClusters_count = (self.Well_DeviatedLength_m*self.Well_ProducerProportion_ratio/self.Stimulation_TargetSpacing_m)-1
        
        #integer datatypes
        self.Well_ProducerCount_wells = int(self.Well_ProducerCount_wells+0.001)
        self.Stimulation_TargetClusters_count = int(self.Stimulation_TargetClusters_count+0.001)
        self.Sales_ContractTime_yr = int(self.Sales_ContractTime_yr+0.001)
        self.Stimulation_Steps_count = int(self.Stimulation_Steps_count+0.001)
        
        #lifespan
        self.Cost_Lifespan_yr = int(0.05+self.Power_LifeSpan_s/yr)
        self.Debt_Tenor_yr = np.min([self.Debt_Tenor_yr,self.Cost_Lifespan_yr-1])
        
        #shear modulus
        self.Rock_ShearModulus_Pa = self.Rock_YoungsModulus_Pa/(2.0*(1.0+self.Rock_PoissonRatio_ratio))
        
        #depth scaling
        if self.Strategy_Target_type.upper() == 'TEMP':
            self.Well_TargetDepth_m = (self.Well_TargetTemp_K-self.Air_Temperature_K)/self.Rock_ThermalGradient_Kpm
        else:
            self.Well_TargetTemp_K = self.Well_TargetDepth_m*self.Rock_ThermalGradient_Kpm+self.Air_Temperature_K #K
        self.Fluid_InjectionDensity_kgpm3 = 1.0/water.v_from_PT(P=self.Circulation_Backpressure_Pa,T=self.Circulation_InjectionTemperature_K)
        self.Fluid_InjectionViscosity_Pas = water.mu_from_T(T=self.Circulation_InjectionTemperature_K)
        self.Fluid_InjectionStatic_Pa = self.Fluid_InjectionDensity_kgpm3*g*self.Well_TargetDepth_m
        self.Fluid_ReservoirDensity_kgpm3 = 1.0/water.v_from_PT(P=self.Fluid_InjectionStatic_Pa,T=self.Well_TargetTemp_K)
        self.Fluid_ReservoirViscosity_Pas = water.mu_from_T(T=self.Well_TargetTemp_K)
        self.Fluid_ReservoirStatic_Pa = 0.5*(self.Fluid_InjectionDensity_kgpm3+self.Fluid_ReservoirDensity_kgpm3)*g*self.Well_TargetDepth_m
        if (self.Fluid_InjectionDensity_kgpm3<500) or (self.Fluid_ReservoirDensity_kgpm3<500):
            print("NOTICE: injected water is not a liquid")
        self.Fluid_HeatCapacity_kJpm3K = ((water.h_from_PT(P=self.Fluid_InjectionStatic_Pa,T=self.Well_TargetTemp_K)
                                          -water.h_from_PT(P=self.Circulation_Backpressure_Pa,T=self.Circulation_InjectionTemperature_K))
                                          /((self.Well_TargetTemp_K-self.Circulation_InjectionTemperature_K)/
                                            (0.5*(self.Fluid_InjectionDensity_kgpm3+self.Fluid_ReservoirDensity_kgpm3))))
        
        #intervals, clusters, and well count (per unit of pattern)
        if self.Strategy_Design_type.upper() == 'AGS':
            self.Well_Intervals_count = 1
            self.Well_Clusters_count = 0
            self.Well_InjectorCount_wells = self.Well_ProducerCount_wells
        elif self.Strategy_Design_type.upper() == 'FGS':
            self.Well_Intervals_count = max([1,int(self.Stimulation_TargetClusters_count)])
            self.Well_Clusters_count = 1
            self.Well_InjectorCount_wells = 1
        elif self.Strategy_Design_type.upper() == 'CGS':
            self.Well_Intervals_count = 1
            self.Well_Clusters_count = max([1,int(self.Stimulation_TargetClusters_count)])
            self.Well_InjectorCount_wells = int(self.Well_Pattern_count*max([1,int(self.Well_ProducerCount_wells-1)]))
        elif self.Strategy_Design_type.upper() in ['O&G','DCM']:
            self.Well_Intervals_count = 1
            self.Well_Clusters_count = max([1,int(self.Stimulation_TargetClusters_count)])
            self.Well_InjectorCount_wells = max([1,int(self.Well_ProducerCount_wells-1)])
        else: # 'EGS'
            self.Well_Intervals_count = 1
            self.Well_Clusters_count = max([1,int(self.Stimulation_TargetClusters_count)])
            self.Well_InjectorCount_wells = 1
        
        #well scaling
        self.Well_ProductionRadius_m = 0.5*self.Well_ProductionDiameter_m
        self.Well_SegmentCount_units = max([1,int(self.Well_SegmentCount_units)])
        
        #stress conditions
        s123 = np.zeros(3,dtype=float)
        s123[0] = self.Rock_Density_kgpm3*g*self.Well_TargetDepth_m
        s123[2] = self.Stress_ShminToSv_ratio*(s123[0]-self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio) + self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio
        s123[1] = self.Stress_ShmaxToShmin_ratio*(s123[2]-self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio) + self.Fluid_ReservoirStatic_Pa*self.Rock_Biot_ratio
        tensor = stress.cauchy()
        tensor.set_sigG_from_Principal(s123[2], s123[1], s123[0], self.Stress_ShminAzimuth_rad, self.Stress_ShminDip_rad)
        self.Tensor_00_Pa = tensor.sigG[0,0]
        self.Tensor_01_Pa = tensor.sigG[0,1]
        self.Tensor_02_Pa = tensor.sigG[0,2]
        self.Tensor_10_Pa = tensor.sigG[1,0]
        self.Tensor_11_Pa = tensor.sigG[1,1]
        self.Tensor_12_Pa = tensor.sigG[1,2]
        self.Tensor_20_Pa = tensor.sigG[2,0]
        self.Tensor_21_Pa = tensor.sigG[2,1]
        self.Tensor_22_Pa = tensor.sigG[2,2]
        str_dip = tensor.get_conjugates(plots=False)
        s123.sort(axis=0)
        self.Stress_S1_Pa = s123[2]
        self.Stress_S2_Pa = s123[1]
        self.Stress_S3_Pa = s123[0]
        self.JointSets_Set1Strike_rad = str_dip[0,0]
        self.JointSets_Set1Dip_rad = str_dip[0,1]
        self.JointSets_Set2Strike_rad = str_dip[1,0]
        self.JointSets_Set2Dip_rad = str_dip[1,1]
        self.JointSets_Set3Strike_rad = str_dip[2,0]
        self.JointSets_Set3Dip_rad = str_dip[2,1]
            
        #stimulation target volume or radius
        self.Stimulation_Steps_count = np.max([1,self.Stimulation_Steps_count])
        self.Stimulation_SeedRadius_m = self.Well_Spacing_m/1.2**(self.Stimulation_Steps_count-1)
        if self.Strategy_Stim_type.upper() in ['RAD','R','RADI','RADIUS']:
            f_R = self.Stimulation_TargetRadius_m #m
            f_w = 8.0*self.Stimulation_CriticalPressure_Pa*(1.0-self.Rock_PoissonRatio_ratio**2.0)*f_R/(np.pi*self.Rock_YoungsModulus_Pa)
            vol = (4.0/3.0)*np.pi*f_R**2.0*0.5*f_w
            self.Stimulation_TargetVolume_m3pfrac = vol*(1.0+self.Stimulation_Leakoff_ratio)
        else: # ['V','VOL','VOLUME']
            f_V = self.Stimulation_TargetVolume_m3pfrac/(1.0+self.Stimulation_Leakoff_ratio)
            vol = np.max([1e-9,f_V])
            f_R = ((3.0/16.0)*vol*self.Rock_YoungsModulus_Pa/(self.Stimulation_CriticalPressure_Pa*(1.0-self.Rock_PoissonRatio_ratio**2.0)))**(1.0/3.0)
            f_w = 8.0*self.Stimulation_CriticalPressure_Pa*(1.0-self.Rock_PoissonRatio_ratio**2.0)*f_R/(np.pi*self.Rock_YoungsModulus_Pa)
            self.Stimulation_TargetRadius_m = f_R
        
        #limit well spacing based on target frac radius
        if self.Strategy_Design_type.upper() in ['EGS','FGS','O&G']:
            self.Well_Spacing_m = np.min([self.Well_Spacing_m, self.Stimulation_TargetRadius_m/1.1])

        #modify stim results using tortuosity for frac efficiency
        if self.Well_Clusters_count > 0:
            self.Stimulation_Volume_m3pfrac = self.Stimulation_TargetVolume_m3pfrac*self.Well_Clusters_count
            self.Well_Clusters_count = max([1,int(self.Well_Clusters_count*(
                self.Fracture_Tortuosity_ratio*(1+stats.norm_trunc(1,0.0,0.5,-1.0,1.0)[0]*self.Strategy_SecondOrderNoise_ratio)))])
            self.Stimulation_Volume_m3pfrac = self.Stimulation_Volume_m3pfrac/self.Well_Clusters_count
        f_V = self.Stimulation_Volume_m3pfrac/(1.0+self.Stimulation_Leakoff_ratio)
        vol = np.max([1e-9,f_V])
        f_R = ((3.0/16.0)*vol*self.Rock_YoungsModulus_Pa/(self.Stimulation_CriticalPressure_Pa*(1.0-self.Rock_PoissonRatio_ratio**2.0)))**(1.0/3.0)
        f_w = 8.0*self.Stimulation_CriticalPressure_Pa*(1.0-self.Rock_PoissonRatio_ratio**2.0)*f_R/(np.pi*self.Rock_YoungsModulus_Pa)
        self.Stimulation_Radius_m = f_R
        self.Stimulation_ApertureCritical_m = f_w
        self.Stimulation_Spacing_m = self.Well_DeviatedLength_m*self.Well_ProducerProportion_ratio/(self.Well_Clusters_count+1)
        
        #simplistic calculations
        #self.Well_DrilledLength_m = self.Well_TargetDepth_m + self.Well_ProducerProportion_ratio*self.Well_DeviatedLength_m*(1 - 0.5*np.sin(self.Well_Dip_rad))
        self.Well_DrilledLength_m = self.Well_TargetDepth_m + self.Well_DeviatedLength_m*(1 - 0.5*np.sin(self.Well_Dip_rad)) 
        self.Stimulation_Pressure_Pa = ((10.7/0.9e-3)*(self.Well_DrilledLength_m*self.Fluid_InjectionViscosity_Pas*self.Fluid_InjectionDensity_kgpm3*g*self.Stimulation_InjectionRate_m3ps**1.852)
                                        /(self.Well_Roughness_metric**1.852*(self.Well_ProductionDiameter_m-2.0*self.Well_WallThickness_m)**4.87) 
                                        +self.Stress_S3_Pa + self.Stimulation_CriticalPressure_Pa - self.Fluid_InjectionStatic_Pa)
        

#**********************************************************************************************************
### save file
#**********************************************************************************************************
def save(s=setup(),timeseries=[],timenames=[],aux=[]):
    out = []
    #inputs
    for var in vars(s):
        if np.isscalar(var):
            out += [[var,getattr(s,var)]]
    #timeseries
    if len(timeseries) > 0:
        for i in range(0,len(timeseries)):
            pkg = '%.3e' %(timeseries[i][0])
            for t in range(1,len(timeseries[i])):
                pkg += ';%.3e' %(timeseries[i][t])
            out += [[timenames[i],pkg]]
    #auxillary
    if len(aux) > 0:
        out += aux
        
    #to file
    csv.save_csv(out=out,filename='setup_payoff.csv',append=True)
    
#**********************************************************************************************************
### convert from GeoDT 5 to GeoDT 6
#**********************************************************************************************************
def convert_5_to_6(filename='*.gts'):
    #load file
    names, data = csv.load_csv(filename)
    
    #match fields
    s = setup(run=False)
    for item in vars(s):
        try:
            #spcial cases
            if item == 'Strategy_DiscreteFractures_path':
                setattr(s,item,['','','','none'])
            
            #match where names are in common
            elif item in names:
                setattr(s,item,data[item])
                
            #port when possible
            elif item == 'Circulation_TargetRate_m3ps':
                o0 = data['Circulation_InjectionRate_m3ps'][0]
                o1 = data['Circulation_InjectionRate_m3ps'][1]
                o2 = data['Circulation_InjectionRate_m3ps'][2]
                o3 = data['Circulation_InjectionRate_m3ps'][3]
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Well_ProductionDiameter_m':
                o = np.copy(data['Well_CasingRadius_m'])
                for i in range(0,3):
                    o[i] = '%s'%(float(o[i])*2)
                setattr(s,item,[o[0],o[1],o[2],o[3]])
            elif item == 'Well_WallThickness_m':
                o = np.copy(data['Well_CasingRadius_m'])
                o2 = np.copy(data['Well_HydraulicRadius_m'])
                for i in range(0,3):
                    o[i] = '%s'%(float(o[i])-float(o2[i]))
                setattr(s,item,[o[0],o[1],o[2],o2[3]])
            elif item == 'Well_AnnularThickness_m':
                o2 = np.copy(data['Well_CasingRadius_m'])
                o = np.copy(data['Well_BoreRadius_m'])
                for i in range(0,3):
                    o[i] = '%s'%(float(o[i])-float(o2[i]))
                setattr(s,item,[o[0],o[1],o[2],o[3]])
            elif item == 'Stimulation_TargetVolume_m3pfrac':
                o = np.copy(data['Stimulation_TargetVolume_m3'])
                o2 = np.copy(data['Stimulation_Clusters_clusters'])
                for i in range(0,3):
                    o[i] = '%s'%(float(o[i])/float(o2[1]))
                setattr(s,item,[o[0],o[1],o[2],o[3]])
            elif item == 'Stimulation_CriticalPressure_Pa':
                o = np.copy(data['Domain_PressureIncrement_Pa'])
                o[1] = '%s'%(float(o[1])*10)
                setattr(s,item,['',o[1],'',o[3]])
            elif item == 'Stress_ShminToSv_ratio':
                o = np.copy(data['StressMin_Coefficient_ratio'])
                setattr(s,item,[o[0],o[1],o[2],o[3]])
            elif item == 'Stress_ShmaxToShmin_ratio':
                o = np.copy(data['StressMin_Coefficient_ratio'])
                o2 = np.copy(data['StressMax_Coefficient_ratio'])
                for i in range(0,3):
                    o[i] = '%s'%(float(o2[i])/float(o[i]))
                setattr(s,item,[o[2],o[1],o[0],o[3]])
            elif item == 'Stress_ShminAzimuth_rad':
                setattr(s,item,data['StressMin_Azimuth_rad'])
            elif item == 'Stress_ShminDip_rad':
                setattr(s,item,data['StressMin_Dip_rad'])
            elif item == 'Stress_Variance_rad':
                setattr(s,item,data['StressMin_Variance_rad'])
            elif item == 'Shear_FrictionCoefficient_ratio':
                o0 = '%s'%(np.tan(float(data['Shear_FrictionAngleMin_rad'][1])))
                o1 = '%s'%(np.tan(float(data['Shear_FrictionAngleNom_rad'][1])))
                o2 = '%s'%(np.tan(float(data['Shear_FrictionAngleMax_rad'][1])))
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Shear_Cohesion_Pa':
                o0 = data['Shear_CohesionMin_Pa'][1]
                o1 = data['Shear_CohesionNom_Pa'][1]
                o2 = data['Shear_CohesionMax_Pa'][1]
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Shear_SlipLength_ratio':
                o0 = data['Fracture_SlipLengthMin_ratio'][1]
                o1 = data['Fracture_SlipLengthNom_ratio'][1]
                o2 = data['Fracture_SlipLengthMax_ratio'][1]
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Shear_DilationSlip_ratio':
                o0 = data['Fracture_DilationSlipMin_ratio'][1]
                o1 = data['Fracture_DilationSlipNom_ratio'][1]
                o2 = data['Fracture_DilationSlipMax_ratio'][1]
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Shear_HydraulicDilation_ratio':
                o0 = data['Fracture_HydraulicDilationMin_ratio'][1]
                o1 = data['Fracture_HydraulicDilationNom_ratio'][1]
                o2 = data['Fracture_HydraulicDilationMax_ratio'][1]
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Tension_FrictionCoefficient_ratio':
                o0 = '%s'%(np.tan(float(data['Tension_FrictionAngle_rad'][0])))
                o1 = '%s'%(np.tan(float(data['Tension_FrictionAngle_rad'][1])))
                o2 = '%s'%(np.tan(float(data['Tension_FrictionAngle_rad'][2])))
                o3 = data['Tension_FrictionAngle_rad'][3]
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Tension_FractureCompressibility_1pPa':
                o0 = data['Fracture_CompressibilityMin_1pPa'][1]
                o1 = data['Fracture_CompressibilityNom_1pPa'][1]
                o2 = data['Fracture_CompressibilityMax_1pPa'][1]
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Stimulation_ProppantCompressibility_1pPa':
                o0 = data['Proppant_CompressibilityMin_1pPa'][1]
                o1 = data['Proppant_CompressibilityNom_1pPa'][1]
                o2 = data['Proppant_CompressibilityMax_1pPa'][1]
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'JointSets_InitialDilation_m':
                o = data['Fracture_InitialHydraulicApertureNom_m']
                o2 = data['Fracture_HydraulicDilationNom_ratio']
                for i in range(0,3):
                    o[i] = '%s'%(float(o[i])/float(o2[1]))
                setattr(s,item,[o[0],o[1],o[2],o[3]])
            elif item == 'Fracture_Tortuosity_ratio':
                o0 = ''
                o1 = data['Fracture_TortuosityMin_ratio'][1]
                o2 = ''
                o3 = 'none'
                setattr(s,item,[o0,o1,o2,o3])
            elif item == 'Stimulation_ProppantPermeability_m2':
                o0 = data['Proppant_InitialPermeabilityMin_m2'][1]
                o1 = data['Proppant_InitialPermeabilityNom_m2'][1]
                o2 = data['Proppant_InitialPermeabilityMax_m2'][1]
                o3 = 'uniform'
                setattr(s,item,[o0,o1,o2,o3])
    
            elif item == 'Sales_ElectricWholesale_USDpkWh':
                setattr(s,item,data['Economics_ElectricitySales_USDpkWh'])
            elif item == 'Sales_ElectricContract_USDpkWh':
                setattr(s,item,data['Economics_ElectricitySales_USDpkWh'])
            elif item == 'Sales_ElectricPurchase_USDpkWh':
                setattr(s,item,data['Economics_ElectricitySales_USDpkWh'])
            
            #use defaults when not found
            else:
                print('not found: %s' %(item))
                setattr(s,item,['',getattr(s,item),'','none'])
        except:
            print('missing input for %s' %(item))
            setattr(s,item,['',getattr(s,item),'','none'])
    
    #save file
    out = []
    for item in vars(s):
        out += [[item, getattr(s,item)[0], getattr(s,item)[1], getattr(s,item)[2], getattr(s,item)[3]]]
    csv.save_set(out,filename[:-4]+'.set')
    
#**********************************************************************************************************
### replace with additions to setup()
#**********************************************************************************************************
class payoff:
    def __init__(self):
        #keymost parameters
        self.Profit_ProjectIRR_ratio = 0.0
        self.Profit_LCOE_USDpMWh = 0.0
        self.Profit_LCOH_USDpMWh = 0.0
        self.Power_RunningNet_kW = 0.0
        self.Power_ProductiveLife_yr = 0.0
        self.Power_RecoveryFactor_ratio = 0.0
        self.Cost_Capital_USD = 0.0
        self.Cost_Operations_USD = 0.0
        self.Land_PowerDensity_Wpm2 = 0.0
        self.Circulation_Recovery_ratio = 0.0
        
        #running power
        self.Power_RunningGross_kW = 0.0
        self.Power_RunningThermal_kW = 0.0
        self.Power_RunningPump_kW = 0.0
        
        # #temperatures
        # self.Temperature_Injection_C = 0.0
        # self.Temperature_Reservoir_C = 0.0
        # self.Temperature_Decline_Cpyr = 0.0
        
        #flow
        # self.Stimulation_Pressure_Pa = 0.0
        # self.Circulation_Pressure_Pa = 0.0
        # self.Rate_Injection_kgps = 0.0
        # self.Rate_Production_kgps = 0.0
        # self.Rate_Leakoff_kgps = 0.0
        
        #economics
        self.Cost_Grants_USD = 0.0
        self.Cost_Credits_USD = 0.0
        self.Cost_Pumping_USD = 0.0
        self.Cost_Drilling_USD = 0.0 
        self.Cost_Stimulation_USD = 0.0 
        self.Cost_Facilities_USD = 0.0 
        self.Cost_Gridconnect_USD = 0.0
        self.Cost_Capital_yr = 0.0 
        self.Cost_Seismic_USD = 0.0
        self.Profit_NetPresentValue_USD = 0.0
        self.Profit_EquityIRR_ratio = 0.0
        self.Profit_SeriesAIRR_ratio = 0.0
        self.Profit_SeriesBIRR_ratio = 0.0
        self.Sales_Electric_USD = 0.0
        self.Sales_Heat_USD = 0.0
        
        #seimsicity
        self.Stress_S3_Pa = 0.0
        self.Stress_S2_Pa = 0.0
        self.Stress_S1_Pa = 0.0
        self.Fluid_ReservoirStatic_Pa = 0.0
        self.Seismicity_MaxQuake_Mw = 0.0
        self.Seismicity_NumQuake_ea = 0.0
        self.Fractures_ShearStim_ea = 0.0
        self.Fractures_HydroStim_ea = 0.0
        self.Fractures_Hydroprop_ea = 0.0
        
        #tidbits
        # self.Drilling_TotalLength_m = 0.0
        # self.Drilling_CementVolume_m3 = 0.0
        # self.Drilling_CasingWeight_kgpm = 0.0
        # self.Stimulation_SandMass_kg = 0.0
        # self.Stimulation_WaterVolume_m3 = 0.0
        # self.Circulation_WorkingVolume_m3 = 0.0
        # self.Circulation_Recovery_pct = 0.0
        # self.Circulation_MakeupWater_m3 = 0.0
        self.CO2_Offset_MtCO2eq = 0.0
        
#**********************************************************************************************************
### testing
#**********************************************************************************************************
if False:
    s = setup()
    save(s)
        
        
"""
        
        
        
        
        
    #     #Emissions offsetting
    #     1e-9*np.sum(self.Power_Thermal_kWh*(self.CO2_Heat_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep 
    #                                 +(self.Power_Gross_kWh+self.Power_Pumping_kWh)*(self.CO2_Electric_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep)
        
                
    #         y = data['Injection_Rate_m3ps']; yname = 'Circulation Rate (m3/s)'
    # # y = data['Rock_ThermalConductivity_WpmK']; yname = 'Thermal Conductivity (W/m-K)'
    # # y = data['Stimulation_Clusters_clusters']; yname = 'Number of Fractures (stages)' ; ylog = False
    
    # # x = data['Well_Spacing_m']; name = 'Well Spacing (m)' ; xlog = False
    # # x = data['Stimulation_ProppantConcentration_m3pm3']; name = 'Proppant Concentration (m3/m3)' ; xlog = False
    # # x = data['Proppant_InitialPermeabilityNom_m2']; name = 'Proppant Permeability (m2)' ; xlog = True
    # # x = data['Stimulation_Clusters_clusters']; name = 'Number of Fractures (stages)' ; xlog = False
    # # x = data['Well_DeviatedLength_m']; name = 'Lateral Length (m)' ; xlog = False
    # # x = data['Rock_ThermalGradient_Kpm']*1e3; name = 'Thermal Gradient (C/km)' ; xlog = False
    # # x = data['Stimulation_InjectionRate_m3ps']; name = 'Stimulation Rate (m3/s)' ; xlog = True
    # # x = data['Well_BoreRadius_m']; name = 'Borehole Radius (m)' ; xlog = False
    # # x = data['StressMin_Coefficient_ratio']; name = 'Stress Condition (0.5-0.8 NF; 0.8-1.2 SS or RF)' ; xlog = False
    # # producers = data['Well_ProducerCount_wells']
    # # injectors = np.max([producers-1,np.ones(len(data))],axis=0)
    # # x = injectors+producers; name = 'Number of Wells (Injectors & Producers)' ; xlog = False
    # # x = data['Well_TargetDepth_m']; name = 'Depth (m)' ; xlog = False #l3
    # # x = data['Stimulation_TargetRadius_m']; name = 'Fracture Radius (m)' ; xlog = False #l3
    # x = data['Injection_BottomPressure_MPa']-data['Production_BottomPressure_MPa']; name = 'Fault Pressure Loss (MPa)'; xlog = True
    # x = data['Fracture_Conductivity_m2m']*3.324e15; name = 'Fracture Conductivity (mD-ft)'; xlog = True
        """
        
        
        
        