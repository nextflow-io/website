import React, { useState, useEffect } from "react";

const countryMap: { [key: string]: string } = {
  ar: "Argentina",
  au: "Australia",
  bd: "Bangladesh",
  be: "Belgium",
  br: "Brazil",
  ca: "Canada",
  ch: "Switzerland",
  cl: "Chile",
  de: "Germany",
  dk: "Denmark",
  es: "Spain",
  et: "Ethiopia",
  fr: "France",
  gb: "United Kingdom",
  gh: "Ghana",
  gr: "Greece",
  id: "Indonesia",
  in: "India",
  it: "Italy",
  ke: "Kenya",
  kr: "South Korea",
  mx: "Mexico",
  ng: "Nigeria",
  no: "Norway",
  nz: "New Zealand",
  rs: "Serbia",
  se: "Sweden",
  tn: "Tunisia",
  tw: "Taiwan",
  us: "United States",
  za: "South Africa",
};

interface AmbassadorFilterProps {
  onFilterChange: (countries: string[]) => void;
}

const AmbassadorFilter: React.FC<AmbassadorFilterProps> = ({ onFilterChange }) => {
  const [selectedCountries, setSelectedCountries] = useState<string[]>([]);
  const [isDropdownOpen, setIsDropdownOpen] = useState<boolean>(false);

  const countries = Object.entries(countryMap)
    .sort(([, a], [, b]) => a.localeCompare(b, "en"))
    .map(([code, name]) => ({ code, name }));

  const handleCountryToggle = (countryCode: string) => {
    const updatedCountries = selectedCountries.includes(countryCode)
      ? selectedCountries.filter(code => code !== countryCode)
      : [...selectedCountries, countryCode];
    
    setSelectedCountries(updatedCountries);
    onFilterChange(updatedCountries);
  };

  const clearAllFilters = () => {
    setSelectedCountries([]);
    onFilterChange([]);
  };

  const removeCountry = (countryCode: string) => {
    const updatedCountries = selectedCountries.filter(code => code !== countryCode);
    setSelectedCountries(updatedCountries);
    onFilterChange(updatedCountries);
  };

  // Close dropdown when clicking outside
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      const target = event.target as HTMLElement;
      if (!target.closest('.country-filter-dropdown')) {
        setIsDropdownOpen(false);
      }
    };

    if (isDropdownOpen) {
      document.addEventListener('click', handleClickOutside);
    }

    return () => {
      document.removeEventListener('click', handleClickOutside);
    };
  }, [isDropdownOpen]);

  return (
    <div className="mb-8 flex flex-col sm:flex-row gap-4 items-center justify-center">
      <label className="">Filter by country:</label>
      <div className="relative country-filter-dropdown">
        <div
          onClick={() => setIsDropdownOpen(!isDropdownOpen)}
          className="px-3 py-2 border border-nextflow-600 hover:border-nextflow-900 bg-white shadow-sm min-w-[300px] max-w-[300px] cursor-pointer"
        >
          <div className="flex items-center min-h-[24px]">
            <div className="flex items-center gap-1 overflow-x-auto flex-1 mr-2">
              {selectedCountries.map((countryCode) => {
                const country = countries.find(c => c.code === countryCode);
                return (
                  <div
                    key={countryCode}
                    className="flex items-center bg-nextflow-100 border border-nextflow-300 rounded-full px-2 py-1 text-xs flex-shrink-0"
                    onClick={(e) => e.stopPropagation()} 
                  >
                    <img
                      className="w-3 h-2 mr-1 rounded overflow-hidden object-cover"
                      src={`https://flagicons.lipis.dev/flags/4x3/${countryCode}.svg`}
                      alt={`${country?.name} flag`}
                    />
                    <span className="mr-1 whitespace-nowrap">{country?.name}</span>
                    <button
                      onClick={(e) => {
                        e.stopPropagation();
                        removeCountry(countryCode);
                      }}
                      className="text-nextflow-600 hover:text-nextflow-900 ml-1"
                      title={`Remove ${country?.name}`}
                    >
                      <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                        <path d="M18 6L6 18M6 6l12 12" />
                      </svg>
                    </button>
                  </div>
                );
              })}
              
              {/* Placeholder text when no selection */}
              {selectedCountries.length === 0 && (
                <span className="text-sm text-gray-500 whitespace-nowrap">Select countries...</span>
              )}
            </div>
            
            {/* Dropdown arrow - always visible on the right */}
            <div className="flex-shrink-0">
              <svg className={`w-4 h-4 transition-transform ${isDropdownOpen ? 'rotate-180' : ''}`} fill="none" stroke="currentColor" viewBox="0 0 24 24">
                <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M19 9l-7 7-7-7" />
              </svg>
            </div>
          </div>
        </div>

        {isDropdownOpen && (
          <div className="absolute top-full left-0 right-0 mt-1 bg-white border border-nextflow-600 shadow-lg min-h-[300px] overflow-y-auto z-10">
            <div className="p-2 border-b">
              {/* <button
                onClick={clearAllFilters}
                className="text-sm text-nextflow-600 hover:text-nextflow-900"
              >
                Clear all filters
              </button> */}
            </div>
            {countries.map(({ code, name }) => (
              <div
                key={code}
                className="flex items-center px-3 py-2 hover:bg-gray-100 cursor-pointer"
                onClick={() => handleCountryToggle(code)}
              >
                <input
                  type="checkbox"
                  checked={selectedCountries.includes(code)}
                    onChange={() => {}} // Handled by the div onClick
                  className="mr-3 accent-nextflow-600"
                />
                <img
                  className="w-6 h-4 mr-3 rounded overflow-hidden object-cover"
                  src={`https://flagicons.lipis.dev/flags/4x3/${code}.svg`}
                  alt={`${name} flag`}
                />
                <span className="text-sm">{name}</span>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
};

export default AmbassadorFilter;
