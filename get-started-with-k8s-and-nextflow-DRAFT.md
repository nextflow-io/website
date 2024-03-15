---
title: Get started with Kubernetes and Nextflow
date: 2019-04-03
type: post
tags: nextflow,kubernetes
status: draft
author: Evan Floden
icon: evan.jpg
---

### Get started with Kubernetes and Nextflow

The following walk-through takes advantage of the recent addition of [kubernetes to docker](https://www.docker.com/kubernetes). As of writing, Kubernetes support in Docker for Mac part of the _edge_ release but is planned to be incorporated in the main releases soon.

The following code runs is confirmed to run on MacOS (High Sierra) running Docker version 18.03.0-ce-rc4-mac57.

1. Open Docker and ensure it is running with Kubernetes
   <img alt='Docker for mac with kubernetes' width='760' src='/img/docker-for-mac-with-kubernetes.png' style='margin:1em auto'/>

2. Install Nextflow if necessary

```
curl -s https://get.nextflow.io | bash
```

3. Create a Kubernetes persistent volume

```
kubectl create -f https://k8s.io/docs/tasks/configure-pod-container/task-pv-volume.yaml
```

4. Create a Kubernetes persistent volume claim

```
kubectl create -f https://k8s.io/docs/tasks/configure-pod-container/task-pv-claim.yaml
```

5. Run any Nextflow pipeline on your new kubernetes clusters

```
nextflow kuberun nextflow-io/rnatoy -v task-pv-claim:/mnt/
```

```

```
