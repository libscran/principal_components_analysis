<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.8">
  <compound kind="file">
    <name>blocked_pca.hpp</name>
    <path>scran/</path>
    <filename>blocked__pca_8hpp.html</filename>
    <class kind="struct">scran::blocked_pca::Options</class>
    <class kind="struct">scran::blocked_pca::Results</class>
    <namespace>scran</namespace>
    <namespace>scran::blocked_pca</namespace>
  </compound>
  <compound kind="file">
    <name>scran.hpp</name>
    <path>scran/</path>
    <filename>scran_8hpp.html</filename>
    <namespace>scran</namespace>
  </compound>
  <compound kind="file">
    <name>simple_pca.hpp</name>
    <path>scran/</path>
    <filename>simple__pca_8hpp.html</filename>
    <class kind="struct">scran::simple_pca::Options</class>
    <class kind="struct">scran::simple_pca::Results</class>
    <namespace>scran</namespace>
    <namespace>scran::simple_pca</namespace>
  </compound>
  <compound kind="struct">
    <name>scran::blocked_pca::Options</name>
    <filename>structscran_1_1blocked__pca_1_1Options.html</filename>
    <member kind="variable">
      <type>bool</type>
      <name>scale</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>a66b6d24de80d8a1b3266ee05921e0774</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>transpose</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>aee11c99fcc8d849abae0b82377bf7a8d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>block_weights::Policy</type>
      <name>block_weight_policy</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>ad8260915fb7563c34e73d37babe12082</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>block_weights::VariableParameters</type>
      <name>variable_block_weight_parameters</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>a2842601c41bb34e500de919b6e321ae5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>components_from_residuals</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>aaa9492f0332ea2e4ddca02e90c6d24fd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>realize_matrix</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>a820f5344683171197213a287e45c1f48</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>ae7cd711769eca1bf813d1089d5309d2c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>irlba::Options</type>
      <name>irlba_options</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Options.html</anchorfile>
      <anchor>a5dbf529efa9137e3a3adf8dd23807bef</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::simple_pca::Options</name>
    <filename>structscran_1_1simple__pca_1_1Options.html</filename>
    <member kind="variable">
      <type>bool</type>
      <name>scale</name>
      <anchorfile>structscran_1_1simple__pca_1_1Options.html</anchorfile>
      <anchor>a7a041f122dfda8e5cd216349f7ffbc86</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>transpose</name>
      <anchorfile>structscran_1_1simple__pca_1_1Options.html</anchorfile>
      <anchor>a92f06efd78601fa38ea5b4fd316b13d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structscran_1_1simple__pca_1_1Options.html</anchorfile>
      <anchor>afc967112d204eac3fb293a79cba93344</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>realize_matrix</name>
      <anchorfile>structscran_1_1simple__pca_1_1Options.html</anchorfile>
      <anchor>aee5d49674f0d28a0abff6e5ee7c04285</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>irlba::Options</type>
      <name>irlba_options</name>
      <anchorfile>structscran_1_1simple__pca_1_1Options.html</anchorfile>
      <anchor>a4b9bc63a775f6f5135c52c341d2efefb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::blocked_pca::Results</name>
    <filename>structscran_1_1blocked__pca_1_1Results.html</filename>
    <templarg>typename EigenMatrix_</templarg>
    <templarg>typename EigenVector_</templarg>
    <member kind="variable">
      <type>EigenMatrix_</type>
      <name>components</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Results.html</anchorfile>
      <anchor>a4652959f3a36f06e538e2ee2ab6164b2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenVector_</type>
      <name>variance_explained</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Results.html</anchorfile>
      <anchor>aedda7568ccba471fc61295e05279e677</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenVector_::Scalar</type>
      <name>total_variance</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Results.html</anchorfile>
      <anchor>a61465d1bd44210464c17f6e5d46e0301</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenMatrix_</type>
      <name>rotation</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Results.html</anchorfile>
      <anchor>a418a2cad1abe790a891c2e4d14a4abfb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenMatrix_</type>
      <name>center</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Results.html</anchorfile>
      <anchor>a30b049698f94e1da057137cfc305d90a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenVector_</type>
      <name>scale</name>
      <anchorfile>structscran_1_1blocked__pca_1_1Results.html</anchorfile>
      <anchor>a690bf5061148806950baec506cf5abc8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>scran::simple_pca::Results</name>
    <filename>structscran_1_1simple__pca_1_1Results.html</filename>
    <templarg>typename EigenMatrix_</templarg>
    <templarg>typename EigenVector_</templarg>
    <member kind="variable">
      <type>EigenMatrix_</type>
      <name>components</name>
      <anchorfile>structscran_1_1simple__pca_1_1Results.html</anchorfile>
      <anchor>a00418c16033d1c6a79c1baecfbe839c4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenVector_</type>
      <name>variance_explained</name>
      <anchorfile>structscran_1_1simple__pca_1_1Results.html</anchorfile>
      <anchor>a5e82333a302bad8b0686653c77af41ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenVector_::Scalar</type>
      <name>total_variance</name>
      <anchorfile>structscran_1_1simple__pca_1_1Results.html</anchorfile>
      <anchor>ad0561bf5b94cb9c081e213fd07d39017</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenMatrix_</type>
      <name>rotation</name>
      <anchorfile>structscran_1_1simple__pca_1_1Results.html</anchorfile>
      <anchor>a6c0d7ecf584567d855343cc2f3c57037</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenVector_</type>
      <name>center</name>
      <anchorfile>structscran_1_1simple__pca_1_1Results.html</anchorfile>
      <anchor>a9714a5458b4242513bc8c1927c2c5d4c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>EigenVector_</type>
      <name>scale</name>
      <anchorfile>structscran_1_1simple__pca_1_1Results.html</anchorfile>
      <anchor>aa783e9c502c32c123c6009f1b394510a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran</name>
    <filename>namespacescran.html</filename>
    <namespace>scran::blocked_pca</namespace>
    <namespace>scran::simple_pca</namespace>
  </compound>
  <compound kind="namespace">
    <name>scran::blocked_pca</name>
    <filename>namespacescran_1_1blocked__pca.html</filename>
    <class kind="struct">scran::blocked_pca::Options</class>
    <class kind="struct">scran::blocked_pca::Results</class>
    <member kind="function">
      <type>Results&lt; EigenMatrix_, EigenVector_ &gt;</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1blocked__pca.html</anchorfile>
      <anchor>a8fd572966cd1661fa93b88b5a05813a0</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; *mat, const Block_ *block, int rank, const Options &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>scran::simple_pca</name>
    <filename>namespacescran_1_1simple__pca.html</filename>
    <class kind="struct">scran::simple_pca::Options</class>
    <class kind="struct">scran::simple_pca::Results</class>
    <member kind="function">
      <type>Results&lt; EigenMatrix_, EigenVector_ &gt;</type>
      <name>compute</name>
      <anchorfile>namespacescran_1_1simple__pca.html</anchorfile>
      <anchor>a4095549335e8e09ad82e91aeb8d44e3c</anchor>
      <arglist>(const tatami::Matrix&lt; Value_, Index_ &gt; *mat, int rank, const Options &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Principal components analysis, duh</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
