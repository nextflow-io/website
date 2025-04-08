import React from "react";
import Button from "../Button";
import Star from "../../../public/img/assets/star.svg";
import nfcore from "../../../public/img/nf-core.png";
import PipelinesIcon from "../../../public/img/assets/pipelinesIcon.svg";
import ContainersIcon from "../../../public/img/assets/containersIcon.svg";

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
      url: "https://seqera.io/pipelines/rnaseq--nf-core/",
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
      url: "https://seqera.io/pipelines/sarek--nf-core/",
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
      url: "https://seqera.io/containers/?packages=bioconda::bcftools=1.2",
    },
    {
      name: "bioconda::samtools",
      description: "Tools for manipulating next-generation sequencing data",
      platforms: ["linux/amd64", "linux/arm64"],
      stars: 380,
      downloads: "95k",
      updated: "3 months ago",
      url: "https://seqera.io/containers/?packages=pip:numpy==2.0.0rc1",
    },
  ];

  return (
    <div className=" mx-auto py-0 px-4">
      <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
        {/* Pipelines Column */}
        <div className="flex flex-col h-full">
          <h2 className="text-lg mb-6 md:text-center flex items-center justify-center">
            <img src={PipelinesIcon.src} alt="Pipelines Icon" className="w-5 h-5 mr-2" />
            Pipelines
          </h2>
          <div className="grid grid-cols-1 auto-rows-fr gap-4">
            {pipelines.map((pipeline, index) => (
              <a 
                key={index} 
                href={pipeline.url} 
                className="border-brand-opacity border p-6 bg-white h-full min-h-[212px] block no-underline transition-all duration-300 hover:border-brand"
              >
                <div className="flex flex-col h-full">
                  <div className="flex items-center gap-2">
                    <img src={nfcore.src} alt="nf-core" className="w-5 h-5" />
                    <span className="text-blue-500 font-medium text-lg">{pipeline.name}</span>
                  </div>

                  <p className=" mt-3 mb-auto text-[14px]">{pipeline.description}</p>

                  <div className="mt-4">
                    <div className="flex flex-wrap gap-2 mb-3">
                      {pipeline.tags.map((tag, tagIndex) => (
                        <span key={tagIndex} className="bg-blue-100 text-blue-700 px-3 py-1 rounded-full text-xs">
                          {tag}
                        </span>
                      ))}
                      {pipeline.additionalTags > 0 && (
                        <span className="text-gray-900 text-sm flex items-center">+ {pipeline.additionalTags} more</span>
                      )}
                    </div>

                    <div className="flex items-center text-gray-900 text-sm gap-3">
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
              </a>
            ))}
          </div>
          <div className="mt-6 flex justify-center">
            <Button url="https://seqera.io/pipelines/" variant="link">
              Launch pipelines
            </Button>
          </div>
        </div>

        {/* Containers Column */}
        <div className="flex flex-col h-full justify-between">
          <h2 className="text-lg mb-6 md:text-center flex items-center justify-center">
            <img src={ContainersIcon.src} alt="Containers Icon" className="w-5 h-5 mr-2" />
            Containers
          </h2>
          <div className="grid grid-cols-1 auto-rows-fr gap-4 h-full">
            {containers.map((container, index) => (
              <a 
                key={index} 
                href={container.url}
                className="border-brand-opacity border p-6 bg-white h-full min-h-[212px] block no-underline transition-all duration-300 hover:border-brand"
              >
                <div className="flex flex-col h-full">
                  <div className="flex items-center gap-2">
                    <span className="text-blue-500 font-medium text-lg flex items-center">
                      <span className="mr-2">
                        <img src={ContainersIcon.src} alt="nf-core" className="w-5 h-5" />
                      </span>
                      {container.name}
                    </span>
                  </div>

                  <p className=" mt-3 mb-auto text-[14px]">{container.description}</p>

                  <div className="mt-4">
                    <div className="flex items-center text-gray-900 text-sm mb-3">
                      {container.platforms.map((platform, platformIndex) => (
                        <span key={platformIndex} className="mr-3">{platform}</span>
                      ))}
                    </div>

                    <div className="flex items-center text-gray-900 text-sm gap-3">
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
              </a>
            ))}
          </div>
          <div className="mt-6 flex justify-center">
            <Button url="https://seqera.io/containers/" variant="link">
              Build containers
            </Button>
          </div>
        </div>
      </div>
    </div>
  );
}


