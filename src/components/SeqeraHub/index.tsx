import React from "react";
import Button from "../Button";
import Star from "../../../public/img/assets/star.svg";
import nfcore from "../../../public/img/nf-core.png";

export default function SeqeraHub() {
  const pipelines = [
    {
      name: "nf-core/rnaseq",
      description:
        "RNA sequencing analysis pipeline using STAR, RSEM, HISAT2 or Salmon with gene/isoform counts and extensive quality control.",
      tags: ["nextflow", "pipeline", "workflow", "nf-core"],
      additionalTags: 2,
      stars: 786,
      size: "198 KB",
      updated: "11 months ago",
    },
    {
      name: "nf-core/sarek",
      description:
        "Analysis pipeline to detect germline or somatic variants (pre-processing, variant calling and annotation) from WGS / targeted sequencing",
      tags: ["nextflow", "pipeline", "workflow", "nf-core"],
      additionalTags: 16,
      stars: 340,
      size: "470 KB",
      updated: "11 months ago",
    },
  ];

  const containers = [
    {
      name: "bioconda::bcftools",
      description: "BCFtools is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF.",
      platforms: ["linux/amd64", "linux/arm64"],
      stars: 450,
      downloads: "120k",
      updated: "2 months ago",
    },
    {
      name: "bioconda::samtools",
      description: "Tools for manipulating next-generation sequencing data",
      platforms: ["linux/amd64", "linux/arm64"],
      stars: 380,
      downloads: "95k",
      updated: "3 months ago",
    },
  ];

  return (
    <div className="max-w-7xl mx-auto py-8 px-4">
      <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
        {/* Pipelines Column */}
        <div>
          <h2 className="text-2xl font-bold mb-6 md:text-center">Pipelines</h2>
          <div className="grid grid-cols-1 auto-rows-fr gap-4">
            {pipelines.map((pipeline, index) => (
              <div key={index} className="border p-6 bg-white shadow-sm h-full">
                <div className="flex flex-col h-full">
                  <div className="flex items-center gap-2">
                    <h3 className="text-blue-500 font-medium text-lg">{pipeline.name}</h3>
                  </div>

                  <p className="text-gray-800 mt-3 mb-auto">{pipeline.description}</p>

                  <div className="mt-4">
                    <div className="flex flex-wrap gap-2 mb-3">
                      {pipeline.tags.map((tag, tagIndex) => (
                        <span key={tagIndex} className="bg-blue-100 text-blue-700 px-3 py-1 rounded-full text-sm">
                          {tag}
                        </span>
                      ))}
                      {pipeline.additionalTags > 0 && (
                        <span className="text-gray-500 text-sm flex items-center">+ {pipeline.additionalTags} more</span>
                      )}
                    </div>

                    <div className="flex items-center text-gray-500 text-sm gap-3">
                      <div className="flex items-center gap-1">
                        <img src={Star.src} alt="Star icon" className="w-4 h-4" />
                        <span>{pipeline.stars}</span>
                      </div>
                      <span>•</span>
                      <span>{pipeline.size}</span>
                      <span>•</span>
                      <span>Updated {pipeline.updated}</span>
                    </div>
                  </div>
                </div>
              </div>
            ))}
          </div>
          <div className="mt-6 flex justify-center">
            <Button url="https://seqera.io/pipelines/" variant="link">
              Browse all pipelines
            </Button>
          </div>
        </div>

        {/* Containers Column */}
        <div>
          <h2 className="text-2xl font-bold mb-6 md:text-center">Containers</h2>
          <div className="grid grid-cols-1 auto-rows-fr gap-4">
            {containers.map((container, index) => (
              <div key={index} className="border p-6 bg-white shadow-sm h-full">
                <div className="flex flex-col h-full">
                  <div className="flex items-center gap-2">
                    <h3 className="text-blue-500 font-medium text-lg flex items-center">
                      <span className="mr-2">
                        <svg xmlns="http://www.w3.org/2000/svg" className="h-5 w-5 text-orange-500" viewBox="0 0 20 20" fill="currentColor">
                          <path fillRule="evenodd" d="M6 5a1 1 0 011-1h6a1 1 0 110 2H7a1 1 0 01-1-1zm-2 3a1 1 0 011-1h10a1 1 0 110 2H5a1 1 0 01-1-1zm0 3a1 1 0 011-1h10a1 1 0 110 2H5a1 1 0 01-1-1zm0 3a1 1 0 011-1h10a1 1 0 110 2H5a1 1 0 01-1-1z" clipRule="evenodd" />
                        </svg>
                      </span>
                      {container.name}
                    </h3>
                  </div>

                  <p className="text-gray-800 mt-3 mb-auto">{container.description}</p>

                  <div className="mt-4">
                    <div className="flex items-center text-gray-500 text-sm mb-3">
                      {container.platforms.map((platform, platformIndex) => (
                        <span key={platformIndex} className="mr-3">{platform}</span>
                      ))}
                    </div>

                    <div className="flex items-center text-gray-500 text-sm gap-3">
                      <div className="flex items-center gap-1">
                        <img src={Star.src} alt="Star icon" className="w-4 h-4" />
                        <span>{container.stars}</span>
                      </div>
                      <span>•</span>
                      <span>{container.downloads} downloads</span>
                      <span>•</span>
                      <span>Updated {container.updated}</span>
                    </div>
                  </div>
                </div>
              </div>
            ))}
          </div>
          <div className="mt-6 flex justify-center">
            <Button url="https://seqera.io/containers/" variant="link">
              Browse all containers
            </Button>
          </div>
        </div>
      </div>
    </div>
  );
}


