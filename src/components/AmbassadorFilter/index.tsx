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
  onFilterChange: (country: string) => void;
}

const AmbassadorFilter: React.FC<AmbassadorFilterProps> = ({ onFilterChange }) => {
  const [selectedCountry, setSelectedCountry] = useState<string>("");

  const countries = Object.entries(countryMap)
    .sort(([, a], [, b]) => a.localeCompare(b, "en"))
    .map(([code, name]) => ({ code, name }));

  const handleCountryChange = (event: React.ChangeEvent<HTMLSelectElement>) => {
    const country = event.target.value;
    setSelectedCountry(country);
    onFilterChange(country);
  };

  return (
    <div className="mb-8 flex flex-col sm:flex-row gap-4 items-center justify-center">
      <label htmlFor="country-filter">Filter by country:</label>
      <div className="relative">
        <select
          id="country-filter"
          value={selectedCountry}
          onChange={handleCountryChange}
          className="px-4 py-2 pl-12 border border-nextflow-600 hover:border-nextflow-900 bg-white shadow-sm min-w-[200px] appearance-none"
        >
          <option value="">All countries</option>
          {countries.map(({ code, name }) => (
            <option key={code} value={code}>
              {name}
            </option>
          ))}
        </select>

        {selectedCountry && (
          <div className="absolute left-3 top-1/2 transform -translate-y-1/2 w-6 h-4 rounded overflow-hidden">
            <img
              className="w-full h-full object-cover"
              src={`https://flagicons.lipis.dev/flags/4x3/${selectedCountry}.svg`}
              alt={`${selectedCountry} flag`}
            />
          </div>
        )}

        {!selectedCountry && (
          <div className="absolute left-3 top-7  transform -translate-y-1/2 pointer-events-none">
            <svg className="w-5 h-5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M21 21l-6-6m2-5a7 7 0 11-14 0 7 7 0 0114 0z"
              />
            </svg>
          </div>
        )}
        {selectedCountry && (
          <div className="absolute right-3 top-2 m-auto flex justify-end h-5 rounded-full">
            <button
              onClick={() => {
                setSelectedCountry("");
                onFilterChange("");
              }}
              className="px-1 py-1  text-sm bg-gray-200 hover:bg-gray-300 rounded-full transition-colors flex items-center justify-center"
              title="Clear filter"
            >
              <svg
                width="12"
                height="12"
                viewBox="0 0 24 24"
                fill="none"
                stroke="currentColor"
                strokeWidth="2"
                strokeLinecap="round"
                strokeLinejoin="round"
              >
                <path d="M18 6L6 18M6 6l12 12" />
              </svg>
            </button>
          </div>
        )}
      </div>
    </div>
  );
};

export default AmbassadorFilter;
