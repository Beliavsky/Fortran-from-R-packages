
# Version 1.0.1: June 10, 2025
  * New parameters:
    - `survivalBenefit`: Generalize survival benefit vectors (previously: unit CF 1 at end of contract)
    - `gammaInZillmer`: As a company-specific feature, include gamma costs (but not beta) in the Zillmer premium
  * Improve test case generation: Also generate code to export sample contract to Excel
  * Add feature to round intermediate value using a `RoundingHelper` class (stored in `Parameters$Hooks$Rounding`)
  * Fix: `InsuranceContract$addBlock` now ensures that unique IDs are used for existing rider contracts 
  

# Version 1.0.0: October 27, 2023
  * Renamed package from LifeInsuranceContracts to LifeInsureR

