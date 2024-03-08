---
title: Selecting the right storage architecture for your Nextflow pipelines
date: 2023-05-04
type: post
description: Exploring the pros and cons of various file systems and cloud object stores for your Nextflow pipeline
image: img/selecting-storage.jpg
tags: nextflow
author: Paolo Di Tommaso
icon: paolo.jpg
---
*In this article we present the various storage solutions supported by Nextflow including on-prem and cloud file systems, parallel file systems, and cloud object stores. We also discuss Fusion file system 2.0, a new high-performance file system that can help simplify configuration, improve throughput, and reduce costs in the cloud.*

At one time, selecting a file system for distributed workloads was straightforward. Through the 1990s, the Network File System (NFS), developed by Sun Microsystems in 1984, was pretty much the only game in town. It was part of every UNIX distribution, and it presented a standard [POSIX interface](https://pubs.opengroup.org/onlinepubs/9699919799/nframe.html), meaning that applications could read and write data without modification. Dedicated NFS servers and NAS filers became the norm in most clustered computing environments.

For organizations that outgrew the capabilities of NFS, other POSIX file systems emerged. These included parallel file systems such as [Lustre](https://www.lustre.org/), [PVFS](https://www.anl.gov/mcs/pvfs-parallel-virtual-file-system), [OpenZFS](https://openzfs.org/wiki/Main_Page), [BeeGFS](https://www.beegfs.io/c/), and [IBM Spectrum Scale](https://www.ibm.com/products/storage-scale-system) (formerly GPFS). Parallel file systems can support thousands of compute clients and deliver more than a TB/sec combined throughput, however, they are expensive, and can be complex to deploy and manage. While some parallel file systems work with standard Ethernet, most rely on specialized low-latency fabrics such as Intel® Omni-Path Architecture (OPA) or InfiniBand. Because of this, these file systems are typically found in only the largest HPC data centers.

## Cloud changes everything

With the launch of [Amazon S3](https://aws.amazon.com/s3/) in 2006, new choices began to emerge. Rather than being a traditional file system, S3 is an object store accessible through a web API. S3 abandoned traditional ideas around hierarchical file systems. Instead, it presented a simple programmatic interface and CLI for storing and retrieving binary objects.

Object stores are a good fit for cloud services because they are simple and scalable to multiple petabytes of storage. Rather than relying on central metadata that presents a bottleneck, metadata is stored with each object. All operations are atomic, so there is no need for complex POSIX-style file-locking mechanisms that add complexity to the design. Developers interact with object stores using simple calls like [PutObject](https://docs.aws.amazon.com/AmazonS3/latest/API/API_PutObject.html) (store an object in a bucket in return for a key) and [GetObject](https://docs.aws.amazon.com/AmazonS3/latest/API/API_GetObject.html) (retrieve a binary object, given a key).

This simple approach was ideal for internet-scale applications. It was also much less expensive than traditional file systems. As a result, S3 usage grew rapidly. Similar object stores quickly emerged, including Microsoft [Azure Blob Storage](https://azure.microsoft.com/en-ca/products/storage/blobs/), [Open Stack Swift](https://wiki.openstack.org/wiki/Swift), and [Google Cloud Storage](https://cloud.google.com/storage/), released in 2010.

## Cloud object stores vs. shared file systems

Object stores are attractive because they are reliable, scalable, and cost-effective. They are frequently used to store large amounts of data that are accessed infrequently. Examples include archives, images, raw video footage, or in the case of bioinformatics applications, libraries of biological samples or reference genomes. Object stores provide near-continuous availability by spreading data replicas across cloud availability zones (AZs). AWS claims theoretical data availability of up to 99.999999999% (11 9's) – a level of availability so high that it does not even register on most [downtime calculators](https://availability.sre.xyz/)!

Because they support both near-line and cold storage, object stores are sometimes referred to as "cheap and deep." Based on current [S3 pricing](https://aws.amazon.com/s3/pricing), the going rate for data storage is USD 0.023 per GB for the first 50 TB of data. Users can "pay as they go" — spinning up S3 storage buckets and storing arbitrary amounts of data for as long as they choose. Some high-level differences between object stores and traditional file systems are summarized below.

<div class="table-responsive">
  <table class="table table-striped table-bordered" style="font-size:80%;">
    <thead>
      <tr>
        <td></td>
        <td>
          <b>Cloud object stores</b>
        </td>
        <td>
          <b>Traditional file systems</b>
        </td>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>
          Interface / access protocol
        </td>
        <td>
          HTTP-based API
        </td>
        <td>
          POSIX interface
        </td>
      </tr>
      <tr>
        <td>
          Cost
        </td>
        <td>
          $
        </td>
        <td>
          $$$
        </td>
      </tr>
      <tr>
        <td>
          Scalability / capacity
        </td>
        <td>
          Practically unlimited
        </td>
        <td>
          Limited
        </td>
      </tr>
      <tr>
        <td>
          Reliability / availability
        </td>
        <td>
          Extremely high
        </td>
        <td>
          Varies
        </td>
      </tr>
      <tr>
        <td>
          Performance
        </td>
        <td>
          Typically lower
        </td>
        <td>
          Varies
        </td>
      </tr>
      <tr>
        <td>
          Support for existing application
        </td>
        <td>
          NO
        </td>
        <td>
          YES
        </td>
      </tr>
    </tbody>
  </table>
</div>

The downside of object storage is that the vast majority of applications are written to work with POSIX file systems. As a result, applications seldom interact directly with object stores. A common practice is to copy data from an object store, perform calculations locally on a cluster node, and write results back to the object store for long-term storage.

## Data handling in Nextflow

Unlike older pipeline orchestrators, Nextflow was built with cloud object stores in mind. Depending on the cloud where pipelines run, Nextflow manages cloud credentials and allows users to provide a path to shared data. This can be a shared file system such as `/my-shared-filesystem/data` or a cloud object store e.g. `s3://my-bucket/data/`.

**Nextflow is exceptionally versatile when it comes to data handling, and can support almost any file system or object store.** Internally, Nextflow uses [executors](https://nextflow.io/docs/latest/executor.html) implemented as plug-ins to insulate pipeline code from underlying compute and storage environments. This enables pipelines to run without modification across multiple clouds regardless of the underlying storage technology.

Suppose an S3 bucket is specified as a location for shared data during pipeline execution. In that case, aided by the [nf-amazon](https://github.com/nextflow-io/nextflow/tree/master/plugins/nf-amazon) plug-in, Nextflow transparently copies data from the S3 bucket to a file system on a cloud instance. Containerized applications mount the local file system and read and write data directly. Once processing is complete, Nextflow copies data to the shared bucket to be available for the next task. All of this is completely transparent to the pipeline and applications.  The same plug-in-based approach is used for other cloud object stores such as Azure BLOBs and Google Cloud Storage.

## The Nextflow scratch directive

The idea of staging data from shared repositories to a local disk, as described above, is not new. A common practice with HPC clusters when using NFS file systems is to use local "scratch" storage.

A common problem with shared NFS file systems is that they can be relatively slow — especially when there are multiple clients. File systems introduce latency, have limited IO capacity, and are prone to problems such as “hot spots” and bandwidth limitations when multiple clients read and write files in the same directory.

To avoid bottlenecks, data is often copied from an NFS filer to local scratch storage for processing. Depending on data volumes, users often use fast solid-state drives or [RAM disks](https://www.mvps.net/docs/how-to-mount-the-physical-memory-from-a-linux-system-as-a-partition/) for scratch storage to accelerate processing.

Nextflow automates this data handling pattern with built-in support for a [scratch](https://nextflow.io/docs/latest/process.html?highlight=scratch#scratch) directive that can be enabled or disabled per process. If scratch is enabled, data is automatically copied to a designated local scratch device prior to processing.

When high-performance file systems such as Lustre or Spectrum Scale are available, the question of whether to use scratch storage becomes more complicated. Depending on the file system and interconnect, parallel file systems performance can sometimes exceed that of local disk. In these cases, customers may set scratch to false and perform I/O directly on the parallel file system.

Results will vary depending on the performance of the shared file system, the speed of local scratch storage, and the amount of shared data to be shuttled back and forth. Users will want to experiment to determine whether enabling scratch benefits pipelines performance.

## Multiple storage options for Nextflow users

Storage solutions used with Nextflow can be grouped into five categories as described below:

- Traditional file systems
- Cloud object stores
- Cloud file systems
- High-performance cloud file systems
- Fusion file system v2.0

The optimal choice will depend on your environment and the nature of your applications and compute environments.

**Traditional file systems** — These are file systems typically deployed on-premises that present a POSIX interface. NFS is the most popular choice, but some users may use high-performance parallel file systems. Storage vendors often package their offerings as appliances, making them easier to deploy and manage. Solutions common in on-prem HPC environments include [Network Appliance](https://www.netapp.com/), [Data Direct Networks](https://www.ddn.com/) (DDN), [HPE Cray ClusterStor](https://www.hpe.com/psnow/doc/a00062172enw), and [IBM Storage Scale](https://www.ibm.com/products/storage-scale-system). While customers can deploy self-managed NFS or parallel file systems in the cloud, most don’t bother with this in practice. There are generally better solutions available in the cloud.

**Cloud object stores** — In the cloud, object stores tend to be the most popular solution among Nextflow users. Although object stores don’t present a POSIX interface, they are inexpensive, easy to configure, and scale practically without limit. Depending on performance, access, and retention requirements, customers can purchase different object storage tiers at different price points. Popular cloud object stores include Amazon S3, Azure BLOBs, and Google Cloud Storage. As pipelines execute, the Nextflow executors described above manage data transfers to and from cloud object storage automatically. One drawback is that because of the need to copy data to and from the object store for every process, performance may be lower than a fast shared file system.

**Cloud file systems** — Often, it is desirable to have a shared file NFS system. However, these environments can be tedious to deploy and manage in the cloud. Recognizing this, most cloud providers offer cloud file systems that combine some of the best properties of traditional file systems and object stores. These file systems present a POSIX interface and are accessible via SMB and NFS file-sharing protocols. Like object stores, they are easy to deploy and scalable on demand. Examples include [Amazon EFS](https://aws.amazon.com/efs/), [Azure Files](https://azure.microsoft.com/en-us/products/storage/files/), and [Google Cloud Filestore](https://cloud.google.com/filestore). These file systems are described as "serverless" and "elastic" because there are no servers to manage, and capacity scales automatically.

Comparing price and performance can be tricky because cloud file systems are highly configurable. For example, [Amazon EFS](https://aws.amazon.com/efs/pricing/) is available in [four storage classes](https://docs.aws.amazon.com/efs/latest/ug/storage-classes.html) – Amazon EFS Standard, Amazon EFS Standard-IA, and two One Zone storage classes – Amazon EFS One Zone and Amazon EFS One Zone-IA. Similarly, Azure Files is configurable with [four different redundancy options](https://azure.microsoft.com/en-us/pricing/details/storage/files/), and different billing models apply depending on the offer selected. To provide a comparison, Amazon EFS Standard costs $0.08 /GB-Mo in the US East region, which is ~4x more expensive than Amazon S3.

From the perspective of Nextflow users, using Amazon EFS and similar cloud file systems is the same as using a local NFS system. Nextflow users must ensure that their cloud instances mount the NFS share, so there is slightly more management overhead than using an S3 bucket. Nextflow users and administrators can experiment with the scratch directive governing whether Nextflow stages data in a local scratch area or reads and writes directly to the shared file system.

Cloud file systems suffer from some of the same limitations as on-prem NFS file systems. They often don’t scale efficiently, and performance is limited by network bandwidth. Also, depending on the pipeline, users may need to stage data to the shared file system in advance, often by copying data from an object store used for long term storage.

For [Nextflow Tower](https://cloud.tower.nf/) users, there is a convenient integration with Amazon EFS. Tower Cloud users can have an Amazon EFS instance created for them automatically via Tower Forge, or they can leverage an existing EFS instance in their compute environment. In either case, Tower ensures that the EFS share is available to compute hosts in the AWS Batch environment, reducing configuration requirements.

**Cloud high-performance file systems** — For customers that need high levels of performance in the cloud, Amazon offers Amazon FSx. Amazon FSx comes in different flavors, including NetApp ONTAP, OpenZFS, Windows File Server, and Lustre. In HPC circles, [FSx for Lustre](https://aws.amazon.com/fsx/lustre/) is most popular delivering sub-millisecond latency, up to 1 TB/sec maximum throughput per file system, and millions of IOPs. Some Nextflow users with data bottlenecks use FSx for Lustre, but it is more difficult to configure and manage than Amazon S3.

Like Amazon EFS, FSx for Lustre is a fully-managed, serverless, elastic file system. Amazon FSx for Lustre is configurable, depending on customer requirements. For example, customers with latency-sensitive applications can deploy FSx cluster nodes with SSD drives. Customers concerned with cost and throughput can select standard hard drives (HDD). HDD-based FSx for Lustre clusters can be optionally configured with an SSD-based cache to accelerate performance. Customers also choose between different persistent file system options and a scratch file system option. Another factor to remember is that with parallel file systems, bandwidth scales with capacity. If you deploy a Lustre file system that is too small, you may be disappointed in the performance.

FSx for Lustre persistent file systems ranges from 125 to 1,000 MB/s/TiB at [prices](https://aws.amazon.com/fsx/lustre/pricing/) ranging from **$0.145** to **$0.600** per GB month. Amazon also offers a lower-cost scratch FSx for Lustre file systems (not to be confused with the scratch directive in Nextflow). At this tier, FSx for Lustre does not replicate data across availability zones, so it is suited to short-term data storage. Scratch FSx for Lustre storage delivers **200 MB/s/TiB**, costing **$0.140** per GB month. This is **~75%** more expensive than Amazon EFS (Standard) and **~6x** the cost of standard S3 storage. Persistent FSx for Lustre file systems configured to deliver **1,000 MB/s/TiB** can be up to **~26x** the price of standard S3 object storage!

**Hybrid Cloud file systems** — In addition to the solutions described above, there are other solutions that combine the best of object stores and high-performance parallel file systems. An example is [WekaFS™](https://www.weka.io/) from WEKA. WekaFS is used by several Nextflow users and is deployable on-premises or across your choice cloud platforms. WekaFS is attractive because it provides multi-protocol access to the same data (POSIX, S3, NFS, SMB) while presenting a common namespace between on-prem and cloud resident compute environments.  Weka delivers the performance benefits of a high-performance parallel file system and optionally uses cloud object storage as a backing store for file system data to help reduce costs.

From a Nextflow perspective, WekaFS behaves like any other shared file system. As such, Nextflow and Tower have no specific integration with WEKA. Nextflow users will need to deploy and manage WekaFS themselves making the environment more complex to setup and manage. However, the flexibility and performance provided by a hybrid cloud file system makes this worthwhile for many organizations.

**Fusion file system 2.0** — Fusion file system is a solution developed by [Seqera Labs](https://seqera.io/fusion) that aims to bridge the gap between cloud-native storage and data analysis workflows. The solution implements a thin client that allows pipeline jobs to access object storage using a standard POSIX interface, thus simplifying and speeding up most operations.

The advantage of the Fusion file system is that there is no need to copy data between S3 and local storage. The Fusion file system driver accesses and manipulates files in Amazon S3 directly. You can learn more about the Fusion file system and how it works in the whitepaper [Breakthrough performance and cost-efficiency with the new Fusion file system](https://seqera.io/whitepapers/breakthrough-performance-and-cost-efficiency-with-the-new-fusion-file-system/).

For sites struggling with performance and scalability issues on shared file systems or object storage, the Fusion file system offers several advantages. [Benchmarks conducted](https://seqera.io/whitepapers/breakthrough-performance-and-cost-efficiency-with-the-new-fusion-file-system/) by Seqera Labs have shown that, in some cases, **Fusion can deliver performance on par with Lustre but at a much lower cost.** Fusion is also significantly easier to configure and manage and can result in lower costs for both compute and storage resources.

## Comparing the alternatives

A summary of storage options is presented in the table below:

<div class="table-responsive">
  <table class="table table-striped table-bordered" style="font-size: 80%;">
    <thead>
      <tr>
        <td></td>
        <td>
          <b>Traditional file systems</b>
        </td>
        <td colspan="3">
          <b>Cloud object storage</b>
        </td>
        <td colspan="3">
          <b>Cloud file systems</b>
        </td>
        <td>
          <b>Fusion FS</b>
        </td>
      </tr>
      <tr>
        <td></td>
        <td>
          NFS, Lustre, Spectrum Scale
        </td>
        <td>
          Amazon S3
        </td>
        <td>
          Azure BLOB storage
        </td>
        <td>
          Google Cloud Storage
        </td>
        <td>
          Amazon EFS
        </td>
        <td>
          Amazon FSX for Lustre
        </td>
        <td>
          Azure File
        </td>
        <td>
          Fusion file system 2.0
        </td>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>
          <b>Deployment model</b>
        </td>
        <td>
          Manual
        </td>
        <td>
          Serverless
        </td>
        <td>
          Serverless
        </td>
        <td>
          Serverless
        </td>
        <td>
          Serverless
        </td>
        <td>
          Serverless
        </td>
        <td>
          Serverless
        </td>
        <td>
          Serverless
        </td>
      </tr>
      <tr>
        <td>
          <b>Access model</b>
        </td>
        <td>
          POSIX
        </td>
        <td>
          Object
        </td>
        <td>
          Object
        </td>
        <td>
          Object
        </td>
        <td>
          POSIX
        </td>
        <td>
          POSIX
        </td>
        <td>
          POSIX
        </td>
        <td>
          POSIX
        </td>
      </tr>
      <tr>
        <td>
          <b>Clouds supported</b>
        </td>
        <td>
          On-prem, any cloud
        </td>
        <td>
          AWS only
        </td>
        <td>
          Azure only
        </td>
        <td>
          GCP only
        </td>
        <td>
          AWS only
        </td>
        <td>
          AWS only
        </td>
        <td>
          Azure only
        </td>
        <td>
          AWS, GCP and Azure <sup>1</sup>
        </td>
      </tr>
      <tr>
        <td>
          <b>Requires block storage</b>
        </td>
        <td>
          Yes
        </td>
        <td>
          Optional
        </td>
        <td>
          Optional
        </td>
        <td>
          Optional
        </td>
        <td>
          Optional
        </td>
        <td>
          No
        </td>
        <td>
          Optional
        </td>
        <td>
          No
        </td>
      </tr>
      <tr>
        <td>
          <b>Relative cost</b>
        </td>
        <td>
          $$
        </td>
        <td>
          $
        </td>
        <td>
          $
        </td>
        <td>
          $
        </td>
        <td>
          $$
        </td>
        <td>
          $$$
        </td>
        <td>
          $$
        </td>
        <td>
          $
        </td>
      </tr>
      <tr>
        <td>
          <b>Nextflow plugins</b>
        </td>
        <td>
          -
        </td>
        <td>
          nf-amazon
        </td>
        <td>
          nf-azure
        </td>
        <td>
          nf-google
        </td>
        <td>
          -
        </td>
        <td>
          -
        </td>
        <td>
          -
        </td>
        <td>
          nf-amazon
        </td>
      </tr>
      <tr>
        <td>
          <b>Tower support</b>
        </td>
        <td>
          Yes
        </td>
        <td>
          Yes, existing buckets
        </td>
        <td>
          Yes, existing BLOB container
        </td>
        <td>
          Yes, existing cloud storage bucket
        </td>
        <td>
          Yes, creates EFS instances
        </td>
        <td>
          Yes, creates FSx for Lustre instances
        </td>
        <td>
          File system created manually
        </td>
        <td>
          Yes, fully automated
        </td>
      </tr>
      <tr>
        <td>
          <b>Dependencies</b>
        </td>
        <td>
          Externally configured
        </td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td>
          Wave Amazon S3
        </td>
      </tr>
      <tr>
        <td>
          <b>Cost model</b>
        </td>
        <td>
          Fixed price on-prem, instance+block storage costs
        </td>
        <td>
          GB per month
        </td>
        <td>
          GB per month
        </td>
        <td>
          GB per month
        </td>
        <td>
          Multiple factors
        </td>
        <td>
          Multiple factors
        </td>
        <td>
          Multiple factors
        </td>
        <td>
          GB per month (uses S3)
        </td>
      </tr>
      <tr>
        <td>
          <b>Level of configuration effort (when used with Tower)</b>
        </td>
        <td>
          High
        </td>
        <td>
          Low
        </td>
        <td>
          Low
        </td>
        <td>
          Low
        </td>
        <td>
          Medium (low with Tower)
        </td>
        <td>
          High (easier  with Tower)
        </td>
        <td>
          Medium
        </td>
        <td>
          Low
        </td>
      </tr>
      <tr>
        <td>
          <b>Works best with:</b>
        </td>
        <td>
          Any on-prem cluster manager (LSF, Slurm, etc.)
        </td>
        <td>
          AWS Batch
        </td>
        <td>
          Azure Batch
        </td>
        <td>
          Google Cloud Batch
        </td>
        <td>
          AWS Batch
        </td>
        <td>
          AWS Batch
        </td>
        <td>
          Azure Batch
        </td>
        <td>
          AWS Batch, Amazon EKS, Azure Batch, Google Cloud Batch <sup>1</sup>
        </td>
      </tr>
    </tbody>
  </table>
</div>

## So what’s the bottom line?

The choice or storage solution depends on several factors. Object stores like Amazon S3 are popular because they are convenient and inexpensive. However, depending on data access patterns, and the amount of data to be staged in advance, file systems such as EFS, Azure Files or FSx for Lustre can also be a good alternative.

For many Nextflow users, Fusion file system will be a better option since it offers performance comparable to a high-performance file system at the cost of cloud object storage. Fusion is also dramatically easier to deploy and manage. [Adding Fusion support](https://nextflow.io/docs/latest/fusion.html) is just a matter of adding a few lines to the `nextflow.config` file.

Where workloads run is also an important consideration. For example, on-premises clusters will typically use whatever shared file system is available locally. When operating in the cloud, you can choose whether to use cloud file systems, object stores, high-performance file systems, Fusion FS, or hybrid cloud solutions such as Weka.

Still unsure what storage solution will best meet your needs? Consider joining our community at [nextflow.slack.com](https://nextflow.slack.com/). You can engage with others, post technical questions, and learn more about the pros and cons of the storage solutions described above.
