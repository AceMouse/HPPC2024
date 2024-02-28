<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="description" content="">
    <meta name="author" content="martin.mortensen@adm.ku.dk">
    <link rel="shortcut icon" href="/assets/ico/favicon.png">

    <title>Login</title>

    <!-- Bootstrap core CSS -->
    <link href="/bootstrap/css/bootstrap.css" rel="stylesheet">

    <!-- Custom styles for this template -->
    <link href="/assets/css/signin.css" rel="stylesheet">
    <link href="/assets/css/ku-login.css" rel="stylesheet">
    <link href="/assets/css/bootstrap-sticky-footer.css" rel="stylesheet">

    <!-- HTML5 shim and Respond.js IE8 support of HTML5 elements and media queries -->
    <!--[if lt IE 9]>
      <script src="/assets/js/html5shiv.js"></script>
      <script src="/assets/js/respond.min.js"></script>
    <![endif]-->
  </head>

  <body>

   <div class="main bg">
    <div class="container" id="login">
     <div class="content">



      <form class="form-signin" action="https://openid.ku.dk/processTrustResult" method="POST">
        
        <h3 class="form-signin-heading">KU OpenID</h3>
        <input type="text" class="form-control" placeholder="KU username" autofocus required name="user">
        <input type="password" class="form-control" placeholder="password" name="pwd" required autocomplete="off">
        
        
        <label class="checkbox">
          <input type="checkbox" value="remember-me" name="remember-me" checked> Remember Trust
        </label>
        <div class="btn-group">
          <button class="btn btn btn-success" type="submit" name="allow">Yes (Allow)</button>
          <button class="btn btn btn-danger " type="submit" name="cancel">No (Cancel)</button>
        </div>
        
        <input type="hidden" name="ct" value="6adec3bff82ffd02e5338b77d301e76ebdd87535">
      </form>

     </div> <!-- /content -->
     <div class="page-content">

    <hr/>



  <!-- Trust root has been validated by OpenID 2 mechanism. -->
<div class="alert info-alert">
  <p>The site <tt>https://erda.dk/</tt> has requested verification
  of your OpenID.</p>
  
</div>





     </div> <!-- /page-content -->
    </div> <!-- /container -->
   </div> <!-- /main bg -->

    <!-- jQuery (necessary for Bootstrap's JavaScript plugins) -->
    <script src="/assets/js/jquery.js"></script>
    <!-- Include all compiled plugins (below), or include individual files as needed -->
    <script src="/assets/bootstrap/js/bootstrap.min.js"></script>

  </body>
</html>
