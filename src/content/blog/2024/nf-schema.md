---
title: "nf-schema: the new and improved nf-validation"
date: 2024-05-01
type: post
description: The post describes nf-schema, the former nf-validation nextflow plugin rewritten and improved.
image: /img/blog-2024-05-01--share.png
tags: nextflow,nf-core,ambassador_post,nf-schema
status: published
author: Nicolas Vannieuwkerke
icon: nicolas.jpeg
community_post: true
ambassador_post: true
---

Check out the revamped Nextflow plugin, nf-schema! It's an enhanced version of nf-validation, utilizing JSON schemas to validate parameters and sample sheets. Unlike its predecessor, it supports the latest JSON schema draft and can convert pipeline-generated files. But what's the story behind its development?

<!-- end-archive-description -->

`nf-validation` is a well-known Nextflow plugin that uses JSON schemas to validate parameters and sample sheets. It can also convert sample sheets to channels using a built-in channel factory. On top of that, it can create a nice summary of pipeline parameters and can even be used to generate a help message for the pipeline.

All of this has made the plugin very popular in the Nextflow community, but it wasn’t without its issues. For example, the plugin uses an older version of the JSON schema draft, namely draft `07` while the latest draft is `2020-12`. It also can’t convert any files/sample sheets created by the pipeline itself since the channel factory is only able to access values from pipeline parameters.

But then `nf-schema` came to the rescue! In this plugin we rewrote large parts of the `nf-validation` code, making the plugin way faster and more flexible while adding a lot of requested features. Let’s see what’s been changed in this new and improved version of `nf-validation`.

# What a shiny new JSON schema draft

To quote the official JSON schema website:

> “JSON Schema is the vocabulary that enables JSON data consistency, validity, and interoperability at scale.”

This one sentence does an excellent job of explaining what JSON schema is and why it was such a great fit for `nf-validation` and `nf-schema`. By using these schemas, we can validate pipeline inputs in a way that would otherwise be impossible. The JSON schema drafts define a set of annotations that are used to set some conditions to which the data has to adhere. In our case, this can be used to determine what a parameter or sample sheet value should look like (this can range from what type of data it has to be to a specific pattern that the data has to follow).

The JSON schema draft `07` already has a lot of useful annotations, but it lacked some special annotations that could elevate our validations to the next level. That’s where the JSON schema draft `2020-12` came in. This draft contained a lot more specialized annotations, like dependent requirements of values (if one value is set, another value also has to be set). Although this example was already possible in `nf-validation`, it was poorly implemented and didn’t follow any consensus specified by the JSON schema team.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-05-01-nfschema-img1a.jpg" alt="meme on bright landscape" />
</div>

# Bye-bye Channel Factory, hello Function

One major shortcoming in the `nf-validation` plugin was the lack of the `fromSamplesheet` channel factory to handle files created by the pipeline (or files imported from another pipeline as part of a meta pipeline). That’s why we decided to remove the `fromSamplesheet` channel factory and replace it with a function called `samplesheetToList` that can be deployed in an extremely flexible way. It takes two inputs: the sample sheet to be validated and converted, and the JSON schema used for the conversion. Both inputs can either be a `String` value containing the path to the files or a Nextflow `file` object. By converting the channel factory to a function, we also decoupled the parameter schema from the actual sample sheet conversion. This means all validation and conversion of the sample sheet is now fully done by the `samplesheetToList` function. In `nf-validation`, you could add a relative path to another JSON schema to the parameter schema so that the plugin would validate the file given with that parameter using the supplied JSON schema. It was necessary to also add this for sample sheet inputs as they would not be validated otherwise. Due to the change described earlier, the schema should no longer be given to the sample sheet inputs because they will be validated twice that way. Last, but certainly not least, this function also introduces the possibility of using nested sample sheets. This was probably one of the most requested features and it’s completely possible right now! Mind that this feature only works for YAML and JSON sample sheets since CSV and TSV do not support nesting.

# Configuration sensation

In `nf-validation`, you could configure how the plugin worked by certain parameters (like `validationSchemaIgnoreParams`, which could be used to exempt certain parameters from the validation). These parameters have now been converted to proper configuration options under the `validation` scope. The `validationSchemaIgnoreParams` has even been expanded into two configuration options: `validation.ignoreParams` and `validation.defaultIgnoreParams`. The former is to be used by the pipeline user to exclude certain parameters from validation, while the latter is to be used by the pipeline developer to set which parameters should be ignored by default. The plugin combines both options so users no longer need to supply the defaults alongside their parameters that need to be ignored.

# But, why not stick to nf-validation?

In February we released an earlier version of these changes as `nf-validation` version `2.0.0`. This immediately caused massive issues in quite some nf-core pipelines (I think I set a new record of how many pipelines could be broken by one release). This was due to the fact that a lot of pipelines didn’t pin the `nf-validation` version, so all these pipelines started pulling the newest version of `nf-validation`. The pipelines all started showing errors because this release contained breaking changes. For that reason we decided to remove the version `2.0.0` release until more pipelines pinned their plugin versions.

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-05-01-nfschema-img1b.jpg" alt="meme on bright landscape" />
</div>

Some discussion arose from this and we decided that version `2.0.0` would always cause issues since a lot of older versions of the nf-core pipelines didn’t pin their nf-validation version either, which would mean that all those older versions (that were probably running as production pipelines) would suddenly start breaking. That’s why there seemed to be only one sensible solution: make a new plugin with the breaking changes! And it would also need a new name. We started collecting feedback from the community and got some very nice suggestions. I made a poll with the 5 most popular suggestions and let everyone vote on their preferred options. The last place was tied between `nf-schemavalidator` and `nf-validationutils`, both with 3 votes. In third place was `nf-checker` with 4 votes. The second place belonged to `nf-validation2` with 7 votes. And with 13 votes we had a winner: `nf-schema`!

So, a fork was made of `nf-validation` that we called `nf-schema`. At this point, the only breaking change was the new JSON schema draft, but some other feature requests started pouring in. That’s the reason why the new `samplesheetToList` function and the configuration options were implemented before the first release of `nf-schema` on the 22nd of April 2024.

And to try and mitigate the same issue from ever happening again, we added an automatic warning when the pipeline is being run with an unpinned version of nf-schema:

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-05-01-nfschema-img1c.png" alt="meme on bright landscape" />
</div>

# So, what’s next?

One of the majorly requested features is the support for nested parameters. The version `2.0.0` already was getting pretty big so I decided not to implement any extra features into it. This is, however, one of the first features that I will try to tackle in version `2.1.0`.

Furthermore, I’d also like to improve the functionality of the `exists` keyword to also work for non-conventional paths (like s3 and azure paths).

It’s also a certainty that some weird bugs will pop up over time, those will, of course, also be fixed.

# Useful links

Here are some useful links to get you started on using `nf-schema`:

If you want to easily migrate from nf-validation to `nf-schema`, you can use the migration guide: https://nextflow-io.github.io/nf-schema/latest/migration_guide/

If you are completely new to the plugin I suggest reading through the documentation: https://nextflow-io.github.io/nf-schema/latest/

If you need some examples, look no further: https://github.com/nextflow-io/nf-schema/tree/master/examples

And to conclude this blog post, here are some very wise words from Master Yoda himself:

<div style="margin-top: 2rem; margin-bottom: 2rem;">
    <img src="/img/blog-2024-05-01-nfschema-img1d.jpg" alt="meme on bright landscape" />
</div>
