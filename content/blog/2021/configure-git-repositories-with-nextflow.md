title=Configure Git private repositories with Nextflow
date=2021-10-21
type=post
description=A step-by-step guide to configure Nextflow with Git hosting solutions.
image=img/nextflow-with-git-hosting.jpg
tags=git,scm
status=published
author=Abhinav Sharma
icon=abhinav.jpg
~~~~~~

Git has become the de-facto standard for source-code version control system and has seen increasing adoption across the spectrum of software development. 

Nextflow provides builtin support for Git and most popular Git hosting platforms such 
as GitHub, GitLab and Bitbucket between the others, which streamline managing versions 
and track changes in your pipeline projects and facilitate the collaboration across 
different users. 

In order to access public repositories Nextflow does not require any special configuration, just use the *http* URL of the pipeline project you want to run 
in the run command, for example: 

```
nextflow run https://github.com/nextflow-io/hello
```

However to allow Nextflow to access private repositories you will need to specifiy 
the repository credentials, and the server hostname in the case of self-managed 
Git server installations.

## Configure access to private repositories

This is done through a file name `scm` placed in the `$HOME/.nextflow/` directory, containing the credentials and other details for accessing a particular Git hosting solution. You can refer to the Nextflow documentation for all the [SCM configuration file](https://www.nextflow.io/docs/edge/sharing.html) options.

All of these platforms have their own authentication mechanisms for Git operations which are captured in the `$HOME/.nextflow/scm` file with the following syntax:

```groovy
providers {
	'<provider-name-1>' {
    		user = value
    		password = value
    		...
	}

	'<provider-name-2>' {
    		user = value
    		password = value
    		...
	}
}
```

Note: Make sure to enclose the provider name with `'` if it containes a `-` or a 
blank character. 

As of the 21.09.0-edge release, Nextflow integrates with the following Git providers:

## GitHub

[GitHub](https://github.com) is one of the most well known Git providers and is home to some of the most popular open-source Nextflow pipelines from the [nf-core](https://github.com/nf-core/) community project.

If you wish to use Nextflow code from a **public** repository hosted on GitHub.com, then you don't need to provide credentials  (`user` and `password`) to pull code from the repository. However, if you wish to interact with a private repository or are running into GitHub API rate limits for public repos, then you must provide elevated access to Nextflow by specifying your credentials in the `scm` file.

It is worth noting that [GitHub recently phased out Git password authentication](https://github.blog/2020-12-15-token-authentication-requirements-for-git-operations/#what-you-need-to-do-today) and now requires that users supply a more secure GitHub-generated *Personal Access Token* for authentication. With Nextflow, you can specify your *personal access token* in the `password` field.

```groovy
providers {
	github {
    		user = 'me'
    		password = 'my-personal-access-token'
	}
}
```

To generate a `personal-access-token` for the GitHub platform, follow the instructions provided  [here](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token). Ensure that the token has at a minimum all the permissions in the `repo` scope.

Once you have provided your username and *personal access token*, as shown above, you can test the integration by pulling the repository code.

```
nextflow pull https://github.com/user_name/private_repo
```

## Bitbucket Cloud

[Bitbucket](https://bitbucket.org/) is a publicly accessible Git solution hosted by Atlassian. Please note that if you are using an on-premises Bitbucket installation, you should follow the instructions for  *Bitbucket Server* in the following section.

If your Nextflow code is in a public Bitbucket repository, then you don't need to specify your credentials to pull code from the repository. However, if you wish to interact with a private repository, you need to provide elevated access to Nextflow by specifying your credentials in the `scm` file.

Please note that Bitbucket Cloud requires your `app password` in the `password` field, which is different from your login password.

```groovy
providers {

	bitbucket {
    		user = 'me'
    		password = 'my-app-password'
	}
}
```

To generate an `app password` for the Bitbucket platform, follow the instructions provided [here](https://support.atlassian.com/bitbucket-cloud/docs/app-passwords/). Ensure that the token has at least `Repositories: Read` permission.

Once these settings are saved in `$HOME/.nextflow/scm`, you can test the integration by pulling the repository code.

```
nextflow pull https://bitbucket.org/user_name/private_repo
```

## Bitbucket Server

[Bitbucket Server](https://www.atlassian.com/software/bitbucket/enterprise) is a Git hosting solution from Atlassian which is meant for teams that require a self-managed solution. If Nextflow code resides in an open Bitbucket repository, then you don't need to provide credentials  to pull code from this repository. However, if you wish to interact with a private repository, you need to give elevated access to Nextflow by specifying your credentials in the `scm` file.

For example, if you'd like to call your hosted Bitbucket server as `mybitbucketserver`, then you'll need to add the following snippet in your `~/$HOME/.nextflow/scm` file.

```groovy
providers {

	mybitbucketserver {
    		platform = 'bitbucketserver'
    		server = 'https://your.bitbucket.host.com'
    		user = 'me'
    		password = 'my-password' // OR "my-token"
	}

}
```

To generate a *personal access token* for Bitbucket Server, refer to the [Bitbucket Support documentation](https://confluence.atlassian.com/bitbucketserver/managing-personal-access-tokens-1005339986.html) from Atlassian.

Once the configuration is saved, you can test the integration by pulling code from a private repository and specifying the `mybitbucketserver` Git provider using the `-hub` option.

```
nextflow pull https://your.bitbucket.host.com/user_name/private_repo -hub mybitbucketserver
```

NOTE: It is worth noting that [Atlassian is phasing out the Server offering](https://www.atlassian.com/migration/assess/journey-to-cloud) in favor of cloud product [bitbucket.org](https://bitbucket.org).

## GitLab

[GitLab](https://gitlab.com) is a popular Git provider that offers features covering various aspects of the DevOps cycle.

If you wish to run a  Nextflow pipeline from a public GitLab repository, there is no need to provide credentials  to pull code. However, if you wish to interact with a private repository, then you must give elevated access to Nextflow by specifying your credentials in the `scm` file.

Please note that you need to specify your *personal access token* in the `password` field.

```groovy
providers {
	mygitlab {
    		user = 'me'
    		password = 'my-password' // or 'my-personal-access-token'
    		token = 'my-personal-access-token'
	}
}
```

In addition, you can specify the `server`  fields for your self-hosted instance of GitLab, by default [https://gitlab.com](https://gitlab.com) is assumed as the server.

To generate a `personal-access-token` for the GitLab platform follow the instructions provided [here](https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html). Please ensure that the token has at least `read_repository`, `read_api` permissions.

Once the configuration is saved, you can test the integration by pulling the repository code using the `-hub` option.

```
nextflow pull https://gitlab.com/user_name/private_repo -hub mygitlab
```

## Gitea

[Gitea server](https://gitea.com/) is an open source Git-hosting solution that can be self-hosted. If you have your Nextflow code in an open Gitea repository, there is no need to specify credentials to pull code from this repository. However, if you wish to interact with a private repository, you can give elevated access to Nextflow by specifying your credentials in the `scm` file.

For example, if you'd like to call your hosted Gitea server `mygiteaserver`, then you'll need to add the following snippet in your `~/$HOME/.nextflow/scm` file.

```groovy
providers {

	mygiteaserver {
    		platform = 'gitea'
    		server = 'https://gitea.host.com'
    		user = 'me'
    		password = 'my-password'
	}
}
```

To generate a *personal access token* for your Gitea server, please refer to the [official guide](https://docs.gitea.io/en-us/api-usage/).

Once the configuration is set, you can test the integration by pulling the repository code and specifying `mygiteaserver` as the Git provider using the `-hub` option.

```
nextflow pull https://git.host.com/user_name/private_repo -hub mygiteaserver
```

## Azure Repos

[Azure Repos](https://azure.microsoft.com/en-us/services/devops/repos/) is a part of Microsoft Azure Cloud Suite. Nextflow integrates natively Azure Repos via the usual `~/$HOME/.nextflow/scm` file.

If you'd like to use the `myazure` alias for the `azurerepos` provider, then you'll need to add the following snippet in your `~/$HOME/.nextflow/scm` file.

```groovy
providers {
	myazure {
    		server = 'https://dev.azure.com'
    		platform = 'azurerepos'
    		user = 'me'
    		token = 'my-api-token'
	}
}
```

To generate a *personal access token* for your Azure Repos integration, please refer to the [official guide](https://docs.microsoft.com/en-us/azure/devops/organizations/accounts/use-personal-access-tokens-to-authenticate?view=azure-devops&tabs=preview-page) on Azure.

Once the configuration is set, you can test the integration by pulling the repository code and specifying `myazure` as the Git provider using the `-hub` option.

```
nextflow pull https://dev.azure.com/org_name/DefaultCollection/_git/repo_name -hub myazure
```

## Conclusion

Git is a popular, widely used software system for source code management. The native integration of Nextflow with various Git hosting solutions is an important feature to facilitate reproducible workflows that enable collaborative development and deployment of Nextflow pipelines.

Stay tuned for more integrations as we continue to improve our support for various source code management solutions!
