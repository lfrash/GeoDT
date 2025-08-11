""" 
Class for calculating cashflow for a geothermal development project:
    1. Debt, equity, fees, and taxes
    2. IRR, LCOE, and NPV
    3. Scheduled power and heat delivery
    4. Two phase pricing schedule (introductory & long term)
    5. Attributes include annual cash flows
"""
#imports
import numpy as np

### economics functions
class cashflow():
    def __init__(self): #defaults to "Durham EGS 3P 3I lpf.xlsx - CF GBP"
        #project parameters
        self.Time_Planning_yr = 1 #yr, final investment decision
        self.Time_Construction_yr = 2 #yr, construction time
        self.Time_Lifespan_yr = 35 #yr, total project life
        self.Cost_Planning_USD = -3.5430875e6 #USD, survey & planning costs
        self.Cost_Reservoir_USD = -129.0447423e6 #USD, drilling and stimulation costs
        self.Cost_Gridconnect_USD = -3.846153846e6 #USD, grid connection costs
        self.Cost_Facilities_USD = -92.1345644e6 #USD, power plant cost
        self.Cost_Grants_USD = 0.0e6 #USD, government grants
        self.Cost_Operations_ratio = 0.02 #ratio, opex vs capex
        #sales pricing
        self.Sales_ElectricWholesale_USDpkWh = 54.645e-3 #USD, market rate sales for electricity
        self.Sales_ElectricContract_USDpkWh = 165.590e-3 #USD, contract sales rate for electricity
        self.Sales_ElectricPurchase_USDpkWh = 291.300e-3 #USD, typical purchase price for customers
        self.Sales_HeatWholesale_USDpkWh = 10.000e-3 #USD, market rate sales for heat
        self.Sales_HeatContract_USDpkWh = 10.000e-3 #USD, contract rate sales for heat
        self.Sales_CapacityFactor_ratio = 0.95 #ratio, capacity factor for plant uptime
        self.Sales_ContractTime_yr = 15 #yr, length of contract for heat and electricity
        #equity, debt, and taxes
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
        #performance, noting that length of these series is shorter than Capex array
        self.Power_Gross_kWh = np.asarray([0,0,0,311220.7353,292546.4382,277972.1724,260977.3087,247599.0859,232420.428,219786.6588,
                              208207.0261,196763.7449,185459.8462,174297.3886,162840.5088,151213.9588,139890.924,
                              128853.0469,118084.7095,107572.5474,97305.06136,86414.00885,75732.30733,65867.86571,
                              56509.57654,47513.06953,38802.21546,30333.07111,22078.4974,14020.82846,6148.07003,
                              -1548.207607,-9074.102275,-16434.13404,-23631.72431])*1e3
        self.Power_Thermal_kWh = np.asarray([0,0,0,0,24382.75325,32510.33766,56893.09091,65020.67532,91689.31169,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,100578.8571,
                                  100578.8571,100578.8571,100578.8571,100578.8571])*1e3
        self.Power_Pumping_kWh = np.asarray([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])*1e3
        #CO2 offset, https://www.gov.uk/government/publications/greenhouse-gas-reporting-conversion-factors-2024
        self.CO2_Electric_kgpkWh = 0.22499 #kgCO2/kWhe, CO2 emitted from 'grid' power
        self.CO2_Heat_kgpkWh = 0.18293 #kgCO2/kWht, CO2 emitted from natural gas for heat 
        self.CO2_Geothermal_kgpkWh = 0.00154
        #exchange rate
        self.Currency_Exchange_USDpGBP = 1.3 #USD, currency exchange rate
        #calculated parameters
        self.re_init()
    def re_init(self):
        #refernce tables
        span = int(self.Time_Lifespan_yr)
        self.Equity_WACC_ratio = self.Equity_Contribution_ratio*self.Equity_CostOfCapital_ratio+(
            (1.0-self.Equity_Contribution_ratio)*self.Debt_InterestRate_ratio*(1.0-self.Tax_Corporate_ratio))
        self.Inflation_Table_series = (1.0+self.Tax_Inflation_ratio)**np.linspace(0,span-1,span)
        self.Discount_Table_series = 1.0/((1.0+self.Equity_WACC_ratio)**np.linspace(0,span-1,span))
        #series
        self.Cost_Capex_USD = np.zeros(span)
        self.Cost_Opex_USD = np.zeros(span)
        a = int(self.Time_Planning_yr)
        b = a + int(self.Time_Construction_yr)
        self.Cost_Capex_USD[a:b] += ((self.Cost_Reservoir_USD+self.Cost_Gridconnect_USD+
                             self.Cost_Facilities_USD+self.Cost_Grants_USD)/(b-a))*np.ones(b-a)
        for i in range(b,span):
            self.Cost_Opex_USD[i] = np.sum(self.Cost_Capex_USD[:i])*self.Cost_Operations_ratio
        self.Cost_Capex_USD[0:a] += (self.Cost_Planning_USD/a)*np.asarray(a)
        #capital and operations costs with inflation
        self.Cost_Capex_USD = self.Cost_Capex_USD*self.Inflation_Table_series
        self.Cost_Opex_USD = self.Cost_Opex_USD*self.Inflation_Table_series
        #depreciation table
        self.Tax_Depreciation_USD = np.zeros(span)
        for i in range(0,span):
            Depreciation_StraightLine_yr = np.max([1,int(1/self.Tax_DepreciationRate_ratio)]) 
            capital = np.ones(Depreciation_StraightLine_yr)*self.Cost_Capex_USD[i]/Depreciation_StraightLine_yr
            j = np.min([i+Depreciation_StraightLine_yr,span])
            self.Tax_Depreciation_USD[i:j] += capital[:j-i]
        #revenue from heat and electricity sales
        j = np.min([span,b+self.Sales_ContractTime_yr])
        self.Sales_Electric_USDpkWh = self.Sales_ElectricWholesale_USDpkWh*self.Inflation_Table_series
        self.Sales_Electric_USDpkWh[:j] = self.Sales_ElectricContract_USDpkWh*self.Inflation_Table_series[:j]
        self.Sales_Heat_USDpkWh = self.Sales_HeatWholesale_USDpkWh*self.Inflation_Table_series
        self.Sales_Heat_USDpkWh[:j] = self.Sales_HeatContract_USDpkWh*self.Inflation_Table_series[:j]
        self.Sales_Energy_USD = self.Power_Gross_kWh*self.Sales_Electric_USDpkWh + self.Power_Thermal_kWh*self.Sales_Heat_USDpkWh
        #pumping losses with electric power and production tax credits
        buy = np.argmin(np.asarray([self.Sales_Electric_USDpkWh,self.Sales_ElectricPurchase_USDpkWh*np.ones(span)]),axis=0)
        self.Cost_Electric_USDpkWh = (self.Sales_ElectricPurchase_USDpkWh*buy+self.Sales_Electric_USDpkWh*(1.0-buy))*self.Inflation_Table_series
        self.Tax_ProductionCredit_USD = (self.Power_Gross_kWh+self.Power_Pumping_kWh*(1.0-buy))*self.Tax_ProductionCredit_USDpkWh
        self.Cost_Energy_USD = self.Power_Pumping_kWh*self.Cost_Electric_USDpkWh
        self.Tax_ProductionCredit_USD[self.Tax_ProductionCredit_USD<0] = 0.0
        self.Tax_CapitalCredit_USD = -1.0*self.Cost_Capex_USD*self.Tax_CapitalCredit_ratio
        #cumulative tax credits (direct pay and liability)
        self.Tax_CumulativeCredit_USD = np.zeros(span)
        self.Tax_DirectPay_USD = np.zeros(span)
        self.Tax_CumulativeCredit_USD[0] = self.Tax_ProductionCredit_USD[0]+(self.Tax_CapitalCredit_USD[0])*(1.0-self.Tax_DirectPay_ratio)
        self.Tax_DirectPay_USD[0] = (self.Tax_CapitalCredit_USD[0])*self.Tax_DirectPay_ratio*(1.0-self.Equity_WACC_ratio)**b
        for i in range(0,span):
            self.Tax_CumulativeCredit_USD[i] = self.Tax_ProductionCredit_USD[i]+(self.Tax_CapitalCredit_USD[i])*(1.0-self.Tax_DirectPay_ratio)+self.Tax_CumulativeCredit_USD[i-1]
            self.Tax_DirectPay_USD[i] = (self.Tax_CapitalCredit_USD[i])*self.Tax_DirectPay_ratio*(1.0-self.Equity_WACC_ratio)**np.min([1,b-i])
        #profit before tax with equity, debt, servicing, interest, and fees
        self.Debt_Equity_USD = np.zeros(span)
        self.Debt_CashFlow_USD = np.zeros(span)
        self.Debt_Cumulative_USD = np.zeros(span)
        self.Debt_Interest_USD = np.zeros(span)
        self.Debt_Reserve_USD = np.zeros(span)
        self.Debt_Arrangement_USD = np.zeros(span)
        self.Profit_CashFlow_USD = np.zeros(span)
        self.Tax_LossCarryForward_USD = np.zeros(span)
        K = self.Cost_Capex_USD+self.Cost_Opex_USD+self.Sales_Energy_USD+self.Cost_Energy_USD+self.Tax_DirectPay_USD
        rB = self.Equity_BrokerFee_ratio
        rA = self.Debt_ArrangementFee_ratio
        rE = self.Equity_Contribution_ratio
        rR = self.Debt_ServiceAccount_ratio
        rI = self.Debt_InterestRate_ratio
        T = self.Debt_Tenor_yr
        x1 = rE*(1.0-rB)
        x2 = 1.0-rB*rE
        self.Debt_CashFlow_USD[0] = 0.0
        self.Debt_Cumulative_USD[0] = 0.0
        self.Debt_Equity_USD[0] = K[0]/(1-rB)
        self.Debt_Interest_USD[0] = 0.0
        self.Debt_Reserve_USD[0] = 0.0
        self.Debt_Arrangement_USD[0] = self.Debt_Equity_USD[0]*rB
        self.Profit_CashFlow_USD[0] = K[0]+self.Debt_Interest_USD[0]+self.Debt_Arrangement_USD[0]
        self.Tax_LossCarryForward_USD[0] = self.Profit_CashFlow_USD[0]-self.Cost_Capex_USD[0]+self.Tax_Depreciation_USD[0]
        for i in range(1,span):
            #debt calculations
            Ci = self.Debt_Cumulative_USD[i-1]
            Ri = self.Debt_Reserve_USD[i-1]
            #development phase
            if i < b:
                self.Debt_CashFlow_USD[i] = ((K[i]+Ci*rI+Ci*rI*rR-Ri)*x2 - (K[i]+Ci*rI)*x1)/(rI*x1+rA*x1-rI*x2-rA*x2+x2-rI*rR*x2)
                self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                self.Debt_Equity_USD[i] = (K[i]+self.Debt_Cumulative_USD[i]*rI+self.Debt_CashFlow_USD[i]*rA)*rE/x2
                self.Debt_Interest_USD[i] = self.Debt_Cumulative_USD[i]*rI 
                self.Debt_Reserve_USD[i] = self.Debt_Interest_USD[i]*rR
            #production phase assuming positive cash flow
            else:
                #debt calculations
                self.Debt_Equity_USD[i] = 0.0
                self.Debt_CashFlow_USD[i] = -1.0*np.min(self.Debt_Cumulative_USD[:i])/T
                self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                self.Debt_Interest_USD[i] = self.Debt_Cumulative_USD[i]*rI 
                self.Debt_Reserve_USD[i] = (self.Debt_Interest_USD[i]-self.Debt_CashFlow_USD[i])*rR
                self.Debt_Arrangement_USD[i] = np.min([0.0,self.Debt_Equity_USD[i]])*rB+np.min([0.0,self.Debt_CashFlow_USD[i]])*rA
                self.Profit_CashFlow_USD[i] = K[i]+self.Debt_Interest_USD[i]+self.Debt_Arrangement_USD[i]
                #payoff loan early if beginning to take on losses other than through capital investment
                if (self.Profit_CashFlow_USD[i-1]<0.0) and (self.Debt_Cumulative_USD[i-1]<-1.0) and (self.Cost_Capex_USD[i-1]>-1.0):
                    self.Debt_Equity_USD[i] = self.Debt_Cumulative_USD[i-1] #-self.Debt_Reserve_USD[i-1]
                    self.Debt_CashFlow_USD[i] = -1.0*self.Debt_Equity_USD[i]
                    self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                #debt is paid at end of tenor
                elif self.Debt_Cumulative_USD[i-1]>-1.0:
                    self.Debt_Equity_USD[i] = 0.0
                    self.Debt_CashFlow_USD[i] = 0.0
                    self.Debt_Cumulative_USD[i] = 0.0
                #production phase if negative cash flow (assumes debt not yet paid)
                elif (self.Profit_CashFlow_USD[i] < 0):
                    self.Debt_Equity_USD[i] = 0.0
                    self.Debt_CashFlow_USD[i] = (K[i]-np.min(self.Debt_Cumulative_USD[:i])/T+Ci*rI+Ci*rI*rR)/(1.0-rI-rA-rI*rR)
                    self.Debt_Cumulative_USD[i] = self.Debt_CashFlow_USD[i]+Ci
                self.Debt_Interest_USD[i] = self.Debt_Cumulative_USD[i]*rI
                if self.Debt_Cumulative_USD[i] < -1.0:
                    self.Debt_Reserve_USD[i] = (self.Debt_Interest_USD[i]-self.Debt_CashFlow_USD[i])*rR
                else:
                    self.Debt_Reserve_USD[i] = self.Debt_Interest_USD[i]*rR
            #cash flow
            self.Debt_Arrangement_USD[i] = np.min([0.0,self.Debt_Equity_USD[i]])*rB+np.min([0.0,self.Debt_CashFlow_USD[i]])*rA
            self.Profit_CashFlow_USD[i] = K[i]+self.Debt_Interest_USD[i]+self.Debt_Arrangement_USD[i]
            self.Tax_LossCarryForward_USD[i] = self.Profit_CashFlow_USD[i]-self.Cost_Capex_USD[i]+self.Tax_Depreciation_USD[i]+np.min([0.0,self.Tax_LossCarryForward_USD[i-1]])
            
        #taxes
        self.Tax_Levied_USD = np.zeros(span)
        self.Profit_BeforeTax_USD = self.Profit_CashFlow_USD-self.Cost_Capex_USD
        self.Tax_LossCarryForward_USD[self.Tax_LossCarryForward_USD>0] = 0.0
        self.Tax_Levied_USD[0] = -(self.Profit_BeforeTax_USD[0]+self.Tax_Depreciation_USD[0])*self.Tax_Corporate_ratio
        self.Tax_Levied_USD[1:] = -(self.Profit_BeforeTax_USD[1:]+self.Tax_Depreciation_USD[1:]+self.Tax_LossCarryForward_USD[:-1])*self.Tax_Corporate_ratio
        self.Tax_Levied_USD[self.Tax_Levied_USD>0] = 0.0
        self.Profit_AfterTax_USD = np.zeros(span)
        for i in range(0,span):
            tax = np.min([0.0,self.Tax_Levied_USD[i]+self.Tax_CumulativeCredit_USD[i]])
            self.Tax_CumulativeCredit_USD[i:] = self.Tax_CumulativeCredit_USD[i:]+tax
            self.Tax_CumulativeCredit_USD[self.Tax_CumulativeCredit_USD<0] = 0.0
            self.Tax_Levied_USD[i] = tax
            self.Profit_AfterTax_USD[i] = self.Profit_BeforeTax_USD[i]+self.Tax_Levied_USD[i]
        
        #cashflows
        self.Profit_Equity_USD = np.zeros(span)
        self.Profit_CashFlowTaxed_USD = self.Profit_CashFlow_USD+self.Tax_Levied_USD
        self.Profit_Equity_USD[0] = self.Profit_AfterTax_USD[0]+self.Cost_Capex_USD[0]-self.Debt_CashFlow_USD[0]+self.Debt_Reserve_USD[0]
        self.Profit_Equity_USD[1:] = self.Profit_AfterTax_USD[1:]+self.Cost_Capex_USD[1:]-self.Debt_CashFlow_USD[1:]+self.Debt_Reserve_USD[1:]-self.Debt_Reserve_USD[:-1]
        keep = np.zeros(span)
        keep[self.Debt_Equity_USD<-1.0] = 1
        keep[self.Profit_AfterTax_USD>0.0] = 1
        keep[self.Debt_Cumulative_USD<-0.1] = 1
        self.Profit_Equity_USD[keep<1] = 0.0
        self.Profit_Company_USD = self.Profit_Equity_USD*self.Equity_CompanyShare_ratio
        if self.Equity_SeriesAShare_ratio > 0:
            self.Profit_SeriesA_USD = self.Profit_Equity_USD*(1.0-self.Equity_CompanyShare_ratio)*self.Equity_SeriesAShare_ratio/(self.Equity_SeriesAShare_ratio+self.Equity_SeriesBShare_ratio)
        else:
            self.Profit_SeriesA_USD = np.zeros(span)
        if self.Equity_SeriesBShare_ratio > 0:
            self.Profit_SeriesB_USD = self.Profit_Equity_USD*(1.0-self.Equity_CompanyShare_ratio)*self.Equity_SeriesBShare_ratio/(self.Equity_SeriesAShare_ratio+self.Equity_SeriesBShare_ratio)
        else:
            self.Profit_SeriesB_USD = np.zeros(span)
        plug = np.ones(b)
        plug[0] = 0.0
        self.Profit_Company_USD[:b] = 0.0
        self.Profit_SeriesA_USD[:b] = (1-plug)*self.Debt_Equity_USD[:b]
        self.Profit_SeriesB_USD[:b] = (plug)*self.Debt_Equity_USD[:b]
        
        #summations
        self.keep = keep
        self.Sales_Heat_USD = np.sum(self.Power_Thermal_kWh*self.Sales_Heat_USDpkWh*keep)
        self.Sales_Electric_USD = np.sum(self.Power_Gross_kWh*self.Sales_Electric_USDpkWh*keep)
        self.Cost_Pumping_USD = np.sum(self.Cost_Energy_USD*keep)
        self.Power_ProductiveLife_yr = np.sum(keep)
        self.Power_RunningGross_kW = np.sum(self.Power_Gross_kWh*keep/(365.25*24))/np.sum(keep)
        self.Power_RunningThermal_kW = np.sum(self.Power_Thermal_kWh*keep/(365.25*24))/np.sum(keep)
        self.Power_RunningPump_kW = np.sum(self.Power_Pumping_kWh*keep/(365.25*24))/np.sum(keep)
        self.Power_RunningNet_kW = np.sum((self.Power_Gross_kWh+self.Power_Pumping_kWh)*keep/(365.25*24))/np.sum(keep)
        
        #IRR and NPV
        self.Profit_EquityIRR_ratio = self.irr(self.Profit_Equity_USD)
        self.Profit_SeriesAIRR_ratio = self.irr(self.Profit_SeriesA_USD)
        self.Profit_SeriesBIRR_ratio = self.irr(self.Profit_SeriesB_USD)
        self.Profit_ProjectIRR_ratio = self.irr(self.Profit_CashFlowTaxed_USD[keep>0])
        self.Profit_NetPresentValue_USD = np.sum(self.Profit_CashFlowTaxed_USD[keep>0]/self.Discount_Table_series[keep>0])
        
        #LCOE and LCOH
        self.Profit_LCOE_USDpMWh = -1e3*(np.sum((self.Cost_Capex_USD+self.Cost_Opex_USD)*self.Discount_Table_series*keep)
                                         /np.max([1,np.sum((self.Power_Gross_kWh+self.Power_Pumping_kWh)*self.Discount_Table_series*keep)]))
        self.Profit_LCOH_USDpMWh = -1e3*(np.sum((self.Cost_Capex_USD+self.Cost_Opex_USD)*self.Discount_Table_series*keep)
                                        /np.max([1,np.sum(self.Power_Thermal_kWh*self.Discount_Table_series*keep)]))
        
        #Emissions offsetting
        self.CO2_Offset_MtCO2eq = 1e-9*np.sum(self.Power_Thermal_kWh*(self.CO2_Heat_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep 
                                    +(self.Power_Gross_kWh+self.Power_Pumping_kWh)*(self.CO2_Electric_kgpkWh-self.CO2_Geothermal_kgpkWh)*keep)
        
    def irr(self,cashflows):
        rate = 0.0
        lo = -10.0
        hi = 100.0
        for i in range(0,25):
            npv = (1 + rate)**np.linspace(0,len(cashflows)-1,len(cashflows))
            npv = np.sum(cashflows/npv)
            if npv > 0:
                lo = rate
            else:
                hi = rate
            rate = (lo + hi)/2
            # print(rate)
        return rate

#example implementation
if False:
    #economics solver configuration
    econ = cashflow()
    
    """
    econ.Currency_Exchange_USDpGBP = ...
    econ.Time_Planning_yr = 1 #yr, final investment decision
    econ.Time_Construction_yr = np.max([2,1+int((Drilling_Time_h+Stimulation_Time_h)/(364.75*24))]) #yr, construction time
    econ.Time_Lifespan_yr = int(econ.Time_Construction_yr+econ.Time_Planning_yr+tp[-1]/(364.75*24*60*60)) #yr, total project life
    econ.Cost_Reservoir_USD = -1.0*(Drilling_Cost_USD+Stimulation_Cost_USD)/econ.Currency_Exchange_USDpGBP #USD, drilling and stimulation costs
    econ.Cost_Facilities_USD = -Power_PlantCost_USD/econ.Currency_Exchange_USDpGBP #USD, power plant cost
    econ.Power_Gross_kWh =   np.asarray(gap + list(0.5*(Bulk[1:]+Bulk[:-1])*me[1:]/mt))*(tp[1]-tp[0])/(60*60)
    econ.Power_Thermal_kWh = np.asarray(gap + list(0.5*(heat[1:]+heat[:-1])))*(tp[1]-tp[0])/(60*60)
    econ.Power_Pumping_kWh = np.asarray(gap + list(0.5*(Pump[1:]+Pump[:-1])))*(tp[1]-tp[0])/(60*60)
    """
    
    econ.re_init()