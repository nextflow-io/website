import React from 'react';
import { AmbassadorNav } from './examples';

const PageExample = () => {
  return (
    <div className="container mx-auto px-4 py-8">
      <h1 className="text-3xl font-bold mb-8 lg:text-center">Ambassador Program</h1>
      
      <div className="flex flex-col lg:flex-row lg:gap-8">
        <div className="mb-6 lg:mb-0">
          <AmbassadorNav />
        </div>
        
        <div className="flex-1">
          <section id="nextflow-ambassadors" className="mb-12">
            <h2 className="text-2xl font-bold mb-4">Nextflow ambassadors</h2>
            <p className="mb-4">
              Nextflow ambassadors are members of the Nextflow community who are passionate about the platform
              and want to help others learn and use it effectively. They are volunteers who contribute to the
              community in various ways, including speaking at events, writing blog posts, and helping others
              on the community forums.
            </p>
            <p>
              Our ambassadors come from diverse backgrounds and have different levels of experience with Nextflow,
              but they all share a common goal: to help others succeed with Nextflow.
            </p>
          </section>
          
          <section id="become-ambassador" className="mb-12">
            <h2 className="text-2xl font-bold mb-4">Become an ambassador</h2>
            <p className="mb-4">
              If you're passionate about Nextflow and want to help others learn and use it effectively,
              consider becoming a Nextflow ambassador. As an ambassador, you'll have the opportunity to:
            </p>
            <ul className="list-disc pl-6 mb-4">
              <li>Share your knowledge and experience with others</li>
              <li>Connect with other members of the Nextflow community</li>
              <li>Attend and speak at Nextflow events</li>
              <li>Receive recognition for your contributions to the community</li>
            </ul>
            <p>
              To become an ambassador, fill out the application form below and tell us about your experience
              with Nextflow and why you want to become an ambassador.
            </p>
          </section>
          
          <section id="our-ambassadors" className="mb-12">
            <h2 className="text-2xl font-bold mb-4">Our ambassadors</h2>
            <p className="mb-4">
              Meet our current ambassadors! These individuals have demonstrated a strong commitment to the
              Nextflow community and are actively helping others succeed with the platform.
            </p>
            <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
              <div className="p-4 border rounded-lg">
                <h3 className="font-bold">Jane Doe</h3>
                <p className="text-sm text-gray-600">Bioinformatics Scientist</p>
              </div>
              <div className="p-4 border rounded-lg">
                <h3 className="font-bold">John Smith</h3>
                <p className="text-sm text-gray-600">Software Engineer</p>
              </div>
              <div className="p-4 border rounded-lg">
                <h3 className="font-bold">Maria Rodriguez</h3>
                <p className="text-sm text-gray-600">Data Scientist</p>
              </div>
            </div>
          </section>
          
          <section id="fundings" className="mb-12">
            <h2 className="text-2xl font-bold mb-4">Fundings</h2>
            <p className="mb-4">
              Nextflow ambassadors can apply for funding to support their community activities,
              such as organizing local meetups, workshops, or other events focused on Nextflow.
            </p>
            <p>
              Funding is available on a case-by-case basis and is intended to help cover expenses
              related to venue rental, refreshments, and other costs associated with organizing community events.
            </p>
          </section>
        </div>
      </div>
    </div>
  );
};

export default PageExample; 